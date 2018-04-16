Kseq cpp interface
==================

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
	while (kseq.read() >= 0)
	{
		std::cout << kseq.name << std::endl; // name, etc. are all std::string
		std::cout << kseq.seq << std::endl;
		if (kseq.qual.size())
			std::cout << "+" << std::endl << kseq.qual << std::endl;
	}
	return 0;
}
```
