/**
 *    \file  arg_parser.cpp
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

#include <arg_parser.hpp>

ArgParser::ArgParser( po::options_description& opt_desc ) : 
  _opt_desc(opt_desc), _vm(po::variables_map())
{
}

ArgParser::~ArgParser(){}

int ArgParser::parse_args( int argc, char **argv ){
  po::positional_options_description p;
  p.add("input-file", -1);
  
  po::store( po::command_line_parser( argc, argv ). 
      options(_opt_desc).positional(p).run(), _vm );
  po::notify(_vm);

  if( _vm.count("help") ){
    std::cout << _opt_desc << std::endl;
    exit(1);
  } 
}

const po::variable_value & ArgParser::operator[](const string& vname){

  if( _vm.count( vname ) ){
    return _vm[vname];
  }else{
    throw runtime_error( vname+" not a valid command line argument " ); 
  }

}
