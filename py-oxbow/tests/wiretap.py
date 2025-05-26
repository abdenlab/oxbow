import asyncio
import functools
import inspect
import pprint
from typing import Any, Final, Literal
from unittest import mock

CALLED: Final[str] = "-> "
RETURNED: Final[str] = "<- "
RAISED: Final[str] = "!! "


def swap_quotes(string):
    if "'" not in string:
        return string
    return string.translate(str.maketrans({"'": '"', '"': "'"}))


class Input:
    """
    Represents a combination of arguments and keyword arguments for
    parameterized tests.

    Examples:
        >>> input = Input(1, 2, a=3, b=4)
        >>> print(input)
        1, 2, a=3, b=4
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return ", ".join(
            (
                *(swap_quotes(repr(v)) for v in self.args),
                *(f"{k}={swap_quotes(repr(v))}" for k, v in self.kwargs.items()),
            )
        )


def obj_to_str(obj):
    try:
        return f"{obj.__module__ or ''}.{obj.__class__.__qualname__}.<object>"
    except AttributeError:
        return pprint.pformat(obj)


class Call:
    """
    Represents a call to a func with its arguments and keyword arguments.

    Examples:
        >>> call = Call("foo", 1, 2, a=3, b=4)
        >>> print(call)
        foo(1, 2, a=3, b=4)
    """

    _serializers = {}

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.input = Input(*args, **kwargs)
        self.calls = []
        self.result: Any = None

    @classmethod
    def register_serializer(cls, type):
        def decorator(func):
            cls._serializers[type] = func
            return func

        return decorator

    def serialize(self) -> str:
        return "\n".join(
            (
                f"{CALLED} {self.func.__module__}.{self.func.__qualname__}({self.input})",
                *(
                    f"    {line}"
                    for c in self.calls
                    for line in c.serialize().splitlines()
                ),
                f"{RAISED if isinstance(self.result, BaseException) else RETURNED} {swap_quotes(self.format(self.result))}",
            )
        )

    @classmethod
    def format(cls, obj):
        if isinstance(obj, (list, tuple)):
            string = type(obj)(cls.format(o) for o in obj[:10])
            if len(obj) > 10:
                string += f"<{len(obj) - 10} more items>"
        elif isinstance(obj, dict):
            string = {k: cls.format(v) for k, v in list(obj.items())[:25]}
            if len(obj) > 25:
                string[f"<{len(obj) - 25} more items>"] = "..."
        else:
            for base in inspect.getmro(type(obj)):
                if serializer := cls._serializers.get(base):
                    string = serializer(obj)
                    break
            else:
                string = obj_to_str(obj)
        if isinstance(string, str) and len(lines := string.splitlines()) > 1:
            return ("\n" + " " * (len(RETURNED) + 1)).join(lines)
        else:
            return string


@Call.register_serializer(BaseException)
def serialize_exception(obj):
    return repr(obj)


def generator():
    yield


@Call.register_serializer(type(generator()))
def serialize_generator(obj):
    return f"{obj.__qualname__}.<generator>"


@Call.register_serializer(type(lambda: None))
def serialize_function(obj):
    return f"{obj.__module__}.{obj.__qualname__}"


class PropertyCall(Call):
    def serialize(self) -> str:
        if self.input.args:
            return super().serialize()
        return "\n".join(
            (
                f"{CALLED} {self.func.__module__}.{self.func.__qualname__}",
                *(
                    f"    {line}"
                    for c in self.calls
                    for line in c.serialize().splitlines()
                ),
                f"{RAISED if isinstance(self.result, BaseException) else RETURNED} {self.format(self.result)}",
            )
        )


class GeneratorCall(Call):
    def serialize(self) -> str:
        if self.input.args:
            return super().serialize()
        return "\n".join(
            (
                f"{CALLED} {self.func.__module__}.{self.func.__qualname__}.<generator>",
                *(
                    f"    {line}"
                    for c in self.calls
                    for line in c.serialize().splitlines()
                ),
                f"{RAISED if isinstance(self.result, BaseException) else RETURNED} {self.format(self.result)}",
            )
        )


class PropertyMock(mock.MagicMock):
    def __get__(self, instance, owner):
        if instance:
            return self(instance, mode="get")

    def __set__(self, instance, value):
        if instance:
            return self(instance, value, mode="set")


class WireTap:
    spying = False
    stack = cursor = []
    spies = []

    def __enter__(self):
        self.spying = True
        for spy in self.spies:
            spy.start()
        return self.stack

    def __exit__(self, exc_type, exc_value, traceback):
        self.stack.clear()
        self.cursor = self.stack
        for spy in self.spies:
            spy.stop()
        self.spies.clear()
        self.spying = False

    def spy(self, obj, key):
        if self.spying:
            raise RuntimeError("Cannot create a spy while the wiretap is active")
        target = getattr(obj, key)
        assert (
            inspect.isfunction(target)
            or inspect.ismethod(target)
            or inspect.iscoroutinefunction(target)
            or isinstance(target, property)
        ), f"Spy target is not a method, function, or property: {type(target)}"
        if isinstance(target, property):
            side_effect = target
            new_callable = PropertyMock
            autospec = None
        else:
            side_effect = target
            new_callable = None
            autospec = True

        def wrapper(*args, **kwargs):
            if not self.spying:
                return side_effect(*args, **kwargs)
            if (
                inspect.ismethod(side_effect)
                or "." in side_effect.__qualname__
                or "self" in inspect.signature(side_effect).parameters
            ):
                call = Call(side_effect, *args[1:], **kwargs)
            else:
                call = Call(side_effect, *args, **kwargs)
            self.cursor.append(call)
            previous = self.cursor
            self.cursor = call.calls
            try:
                result = side_effect(*args, **kwargs)
            except BaseException as e:
                call.result = e
                raise e
            else:
                call.result = result
                if inspect.isgenerator(result):
                    return generator_wrapper(side_effect, result)
                return result
            finally:
                self.cursor = previous

        async def async_wrapper(*args, **kwargs):
            assert inspect.iscoroutinefunction(side_effect)
            if not self.spying:
                return await side_effect(*args, **kwargs)
            if (
                inspect.ismethod(side_effect)
                or "." in side_effect.__qualname__
                or "self" in inspect.signature(side_effect).parameters
            ):
                call = Call(side_effect, *args[1:], **kwargs)
            else:
                call = Call(side_effect, *args, **kwargs)
            self.cursor.append(call)
            previous = self.cursor
            self.cursor = call.calls
            try:
                result = await side_effect(*args, **kwargs)
            except BaseException as e:
                call.result = e
                raise e
            else:
                call.result = result
                if inspect.isgenerator(result):
                    return generator_wrapper(side_effect, result)
                return result
            finally:
                self.cursor = previous

        def property_wrapper(*args, mode: Literal["get", "set"], **kwargs):
            if mode == "get":
                se = side_effect.fget
            else:
                se = side_effect.fset
            assert se
            if not self.spying:
                return se(*args, **kwargs)
            call = PropertyCall(se, *args[1:], **kwargs)
            self.cursor.append(call)
            previous = self.cursor
            self.cursor = call.calls
            try:
                result = se(*args, **kwargs)
            except BaseException as e:
                call.result = e
                raise e
            else:
                call.result = result
                if inspect.isgenerator(result):
                    return generator_wrapper(se, result)
                return result
            finally:
                self.cursor = previous

        def generator_wrapper(func, generator):
            try:
                for value in generator:
                    if not self.spying:
                        yield value
                    call = GeneratorCall(func)
                    self.cursor.append(call)
                    previous = self.cursor
                    self.cursor = call.calls
                    call.result = value
                    yield value
                    self.cursor = previous
            except BaseException as e:
                call.result = e
                raise e
            finally:
                self.cursor = previous

        if asyncio.iscoroutinefunction(target):
            wrapped = functools.update_wrapper(async_wrapper, side_effect)
        elif isinstance(target, property):
            wrapped = functools.update_wrapper(property_wrapper, side_effect)
        else:
            wrapped = functools.update_wrapper(wrapper, side_effect)

        spy_obj = mock.patch.object(
            obj,
            key,
            autospec=autospec,
            new_callable=new_callable,
            side_effect=wrapped,
        )
        spy_obj.spy_return = None
        spy_obj.spy_return_list = []
        spy_obj.spy_exception = None
        self.spies.append(spy_obj)
        return spy_obj
