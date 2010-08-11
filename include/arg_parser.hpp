/**
 *    \file  arg_parser.hpp
 *   \brief  
 *
 *  Detailed description starts here.
 *
 *  \author  Rob Patro (RP), rob@cs.umd.edu
 *
 *  \internal
 *    Created:  04/29/2008
 *   Revision:  $Id: doxygen.templates.example,v 1.4 2007/08/02 14:35:24 mehner Exp $
 *   Compiler:  gcc/g++
 *    Company:  University of Maryland, College Park
 *  Copyright:  Copyright (c) 2008, Rob Patro
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 * =====================================================================================
 */

#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

/** Standard Include */
#include <iostream>
#include <stdexcept>

/** Boost Includes */
#include <boost/program_options.hpp>
#include <boost/utility.hpp>

/** Namespace Uses */
using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::runtime_error;

/** Namespace Declarations */
namespace po = boost::program_options;

class ArgParser : boost::noncopyable
{
  public:
    explicit ArgParser( po::options_description& opt_desc );
    ~ArgParser();

    int parse_args( int argc, char **argv );
    const po::variable_value& operator[](const string& vname);
    
  private:
    po::options_description _opt_desc;  
    po::variables_map _vm;
};
#endif /** ARG_PARSER_HPP */
