Kseq c++ re-implementation
==========================

c++ re-implementation of kseq.h from klib using std::istream

Synopsis
========

This is a re-implementation of kseq.h in c++ dialect.

* re-designed in object-oriented fashion.
* reading function renamed from `kseq_read(kseq_t*)` to `kseq_parser.read()`
* the input stream now completely resides in std::istream, using of `FILE*` is no longer needed (this means `gzFile` should be wrapped into stream also)
* data of parsed seq, namely name, comments, seq, qual are `std::string` instead of `kstring_t*`

Usage
=====

```c++
#include <iostream>
#include "kseq_cpp.h"

int main(const int argc, const char *argv[])
{
	kseq_parser kseq("path/to/fasta.or.fastq.file");
	while (kseq.read_seq() >= 0)
	{
		std::cout << kseq.name << std::endl; // name, etc. are all std::string
		std::cout << kseq.seq << std::endl;
		if (kseq.qual.size())
			std::cout << "+" << std::endl << kseq.qual << std::endl;
	}
	return 0;
}
```

Survive/fail cases
==================

This parser is designed to pass/survive or fail at below specific cases:

### 1. Pass: normal fastq/fasta or hybrid fastq/fasta file

```fasta
>fasta
<---- seq data ---->
```

```
@fastq_seq
<---- seq data ---->
+sample_fastq_seq
<---- qual data --->
```

Parser should pass under theses situations as normal.
Should also pass in case both fasta/fastq seqs are found in a single file.


### 2. Pass: multiple line sequence/quality data

```fastq
@fastq
<---- seq data 1 ---->
<---- seq data 2 ---->
+fastq
<---- qual data 1 --->
<---- qual data 2 --->
<---- qual data 3 --->
```

If either/both seq data and quality data comes in multiple lines (number of lines is not necessary to be the same), the parser should parse normally.


### 3. Pass: missing fastq quality header

```fastq
@fastq
<---- seq data ---->
+
<---- qual data --->
```

The quality data is assumed to be immediately following the sequence data of the same fastq seq.
Hence in case the quality header line contains only one '+', the parser should also function normally.
In fact, if the quality header line does include the identifier, it will be ignored anyway.

### 4. Pass: empty lines

```fastq
empty line
@fastq
empty line
<---- seq data 1 ---->
empty line
<---- seq data 2 ---->
empty line
+
empty line
<---- qual data 1 --->
empty line
<---- qual data 2 --->
empty line
```

If any/all of the positions above contains one or more empty lines, it should parse normally.


### 5. Survive: missing quality data

```fastq
@fastq
<---- seq data ---->

@fastq
<---- seq data ---->
+
```

Above two cases contain fastq seqs without quality data.
The parser should survive under these circumanstances but give warnings (i.e. have return codes greater than 0).
This code will be the same as `kseq_parser::status::no_qual` constant.

### Fail: length of seq and quality data not match

```fastq
@fastq
<---- seq data ---->
+
<-- qual data -->
```
If a fastq file does contain non-zero length quality data, however, its length does not match the seq data,
the parser should fail, and return a negative value `kseq_parser::status::fail`

It is designed in this way as the difficulties in parse Sanger format fastq,
where both allowing the seq and quality leading characters (i.e. '@' and '+') to appear in the quality data, and multi-line quality data simultaneously.

To find an easy solution, the parser assumes any quality data, whenever available, will be in the exact same length as the corresponding seq data (line breakers are not counted).
Hence it will extract the exact same number of characters from the stream when parsing the quality data, no matter what actual character it is (only line breaker is excluded).
In this approach, if the lengths do not match, part of the next seq might be extracted before being actually parsed.
