# Why Oxbow?

Most of the tooling for working with data derived from next-generation and emerging DNA sequencing (NGS) assays is tightly coupled to the specialized ways in which NGS data is stored or queried. 

The result? Most of the business logic of genomic software is consumed with complex and redundant input/output. Genomic workflows shuffle data through meandering, labor-intensive, and time-consuming transformations, and inevitably evolve into <a href="https://xkcd.com/2916/#xt=0&yt=70" target="_blank" rel="noopener">Rube Goldberg machines</a>.

<p align="center">
    <a title="Rube Goldberg, Public domain, via Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File:Rube_Goldberg%27s_%22Self-Operating_Napkin%22_(cropped).gif">
        <img width="256" alt="Self-Operating Napkin" src="https://upload.wikimedia.org/wikipedia/commons/a/a9/Rube_Goldberg%27s_%22Self-Operating_Napkin%22_%28cropped%29.gif?20100331201551">
    </a>
</p>

Can we alleviate the friction imposed by the bioinformatic silo? Since most NGS file formats can be modeled as tabular or "relational", we can take advantage of a modern standard for in-memory representation of tabular data called [Arrow](https://arrow.apache.org/). Many popular open-source data libraries today work directly on columnar Arrow representations with zero copying or transformation.

The goal of Oxbow is to efficiently and flexibly **translate streams of records from conventional NGS files into Arrow representations on demand**. The wide adoption and interoperability of Arrow means that Oxbow makes NGS data accessible to a growing array of downstream, general-purpose and interchangeable data technologies, as opposed to bioinformatics toolsuites designed around a specific file format. Furthermore, Oxbow lets you easily apply these technologies to datasets that may be much larger than memory, for example, by using streaming, lazy evaluation, or distributed computing.

The main aspiration of the project is to help bring genomics into the [composable data systems](https://voltrondata.com/codex) fold.