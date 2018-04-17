Kseq c++ re-implementation
==========================

c++ re-implementation of kseq.h from klib using std::istream

Synopsis
========

This is a re-implementation of kseq.h in c++ dialect.

* re-designed in object-oriented fashion.
* reading function renamed from 'kseq_read(kseq_t\*)' to 'kseq_parser.read()'
* the input stream now completely desides in the std::istream side, no need to used FILE* (this means gzFile should be wrapped into stream also)
* name, comments, seq, qual are all std::string instead of kstring_t*

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

Parser should pass under theses situations as normal. Should also pass in case both fasta/fastq seqs are found in a single file.


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


### Survive: missing quality data

```fastq
@fastq
<---- seq data ---->

@fastq
<---- seq data ---->
+
```

Above two cases contains fastq seqs without quality data.
The parser should survive under these circumanstances but gives out warning (i.e. has return codes greater than 0).

### Fail: length of seq and quality data not match

```fastq
@fastq
<---- seq data ---->
+
<-- qual data -->
```
If a fastq file contains non-zero length quality data, however, the length of this quality data is not the same as the seq data,
the parser should fail, and returns a negative value.
It is designed in this way as the difficulties in parse Sanger format fastq,
where the seq and quality leading character (i.e. '@' and '+') will also appear in the quality data.
The parser, however, assumes quality data and seq data have the exact same length,
thus will extract the exact same number of characters from the stream after pass a quality header line, and treat them as quality data
If the lengths not match, the parse of next seq will be problematic.
The only exception is extra line breakers. They are trimmed when counting 'number of characters'.
