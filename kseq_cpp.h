/*==============================================================================
 * this library is a reimplementation of kseq.h in klib for a more seamless
 * interface within c++ environment
 *   -this implementation resides in c++ and is not a simple objecttive-oriented
 *	    wrapper of the original kseq.h
 *   -the original C lib: https://github.com/attractivechaos/klib/
 *      -> full URL: https://github.com/attractivechaos/klib/blob/master/kseq.h
 *==============================================================================
 * this programs is released under the same lisence as the original lib (MIT)
 * a full copy of lisence should be obtained as in LISENCE along with any
 * distribution of this program
 * IMPORTANT DECLARATION:
 *   -this is a free and open source software;
 *   -any one obtained this software will be granted the rights to copy, modify
 *    and redistribute this software
 *   -there is no garanteed function of this software under any circumastances
 *   -the author should not be liable for any losses by using this software
 *    under any circumastances
 * more details can be found in LISENCE
 *==============================================================================
 * author: Guangyu Li <li.gua@husky.neu.edu>
 *   April 2018
 * original C lib author: Attractive Chaos <attractor@live.co.uk>
 *   05 March 2012 (last modification)
 *==============================================================================
 */

#ifndef __KSEQ_CPP_H__
#define __KSEQ_CPP_H__

#include <limits> // std::numeric_limits
#include <string>
#include <fstream>
#include <memory>


#if (__cplusplus < 201103L)
#define _auto_ptr_t std::auto_ptr
#else
#define _auto_ptr_t std::unique_ptr
#endif


namespace kseq_cpp
{

typedef std::string kstring_t;

class kseq_parser
{
public:
	struct status
	{
		static const int good = 0; // upen successful parse
		static const int no_qual = 1; // when no quality found, still OK
		static const int eof = -1; // reach EOF
		static const int bad_file = -255; // when parse failed
	};

public:
	// can be directly accessed and modified
	// name: seq name (a.k.a. seq definition)
	// comments: contents after the first space in def. line (if any)
	// seq: sequence data
	// qual: quality data (only meaningful for fastq format)
	kstring_t name, comments, seq, qual;

	/*==========================================================================
	 * constructor from file name
	 *==========================================================================
	 * in this case, a ifstream will be opened internally
	 *   -new an ifs object, managed by ifs_prt
	 *   -this->is points to the opened ifs
	 *   -read data only through this->is
	 *   -close ifs when finish (in destructor)
	 *   -delete this ifs in (in destructor)
	 */
	explicit kseq_parser(const char *file_name)
	// this constructor also opens a ifstream instance for input file
	{
		ifs_prt.reset(new std::ifstream(file_name));
		is = ifs_prt.get();
	};

	/*==========================================================================
	 * constructor from istream
	 *==========================================================================
	 * if an input stream is already been opened elsewhere, the only thing need
	 * to do is to work through this istream, and no open/close will be done
	 *   -this->is set to the external istream
	 *      -this stream should be completely managed externally
	 *      -will not be closed by destructor
	 *   -ifs_prt kept as NULL
	 *   -read data only through this->is
	 */
	explicit kseq_parser(std::istream *_stream)
	{
		ifs_prt.reset();
		is = _stream;
	};

	/*==========================================================================
	 * destructor
	 *==========================================================================
	 * destructor deals with the newed ifstream (if did)
	 *   -ifs_prt->close() close stream
	 *   -release ifs_prt
	 */
	~kseq_parser()
	{
		if (ifs_prt.get())
		{
			ifs_prt->close();
			ifs_prt.reset();
		}
		this->is = NULL;
	};

	/*==========================================================================
	 * parse next seq from this->is
	 *==========================================================================
	 * return values:
	 * when >= 0, meaningful seq is parsed
	 *   -status::good (0): succesfully parsed a seq (despite reaching eof)
	 *   -status::no_qual (1): succesfully parsed seq but no quality data found
	 * when < 0, something bad happened and parse failed
	 *   -status::eof (-1): normal exit; any subsequent call after reaching eof
	 *   -status::bad_file (-255): bad file format, parse failed
	 */
	int read_seq(void)
	{
		static char c;
		static bool is_fastq;
		static std::size_t brk_pos;
		static kstring_t line_buf;

		if (is->good())
		{
			while ((is->get(c) && (c != '>') && (c != '@')));
			// decide if is fastq
			is_fastq = (c == '@') ? true : false;

			std::getline(*is, name); // read name from stream
			// for the name don't use line_buf
			if (brk_pos = name.find(' ') != std::string::npos)
			// split comments if any, otherwise empty comments
			{
				comments = name.substr(brk_pos + 1, name.size());
				name.erase(brk_pos);
			}
			else
				comments.clear();
			// read seq
			seq.clear();
			qual.clear();
			while (true)
			{
				if (((c = is->peek()) == EOF) || (c == '>') || (c == '@'))
					// end of current seq
					// if is_fastq, obviously no quality data is found
					// otherwise is good (fasta)
					if (is_fastq)
						return kseq_parser::status::no_qual;
					else
						return kseq_parser::status::good;
				else if (c == '+')
					// proceed to read quality data
					break;
				else
					std::getline(*is, line_buf);
					if (line_buf.size())
						seq += line_buf;
					// otherwise continue read seq data
					continue;
			}
			// here start reading quality data
			// first, ignore this entire line to discard quality def line
			is->ignore(std::numeric_limits< std::streamsize >::max(),
				'\n');
			qual.reserve(seq.size());
			while (qual.size() < seq.size())
			{
				std::getline(*is, line_buf);
				if (line_buf.size())
					qual += line_buf;
			}
			if (qual.size() != seq.size()) // length of seq and qual not match
				return kseq_parser::status::bad_file; // oops! corrupted file
			// discard following '\n'
			while ((c = is->peek()) == '\n')
				is->ignore(1);
			if (((c = is->peek()) == EOF) || (c == '>') || (c == '@'))
				return kseq_parser::status::good; // still good; not eof
			else
				return kseq_parser::status::bad_file; // oops! corrupted file
		}
		else
			return kseq_parser::status::eof;
	};

private:
	_auto_ptr_t< std::ifstream > ifs_prt; // used for pointer managing
	// ifs_prt remains empty if the ifs is created externally

	std::istream *is; // this is used for actually loading data
	kstring_t line_buf; // buffer for reading a single line
};

}; /* namespace kseq_cpp */


#endif /* __KSEQ_CPP_H__ */
