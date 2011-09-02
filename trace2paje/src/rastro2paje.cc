// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Sex, 02 Set 2011 14:34:41 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  rastro2paje.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26-07-2011 13:01:06 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <fstream>
#include <iostream>
#include <algorithm>
#include "semantics.hh"
#include "container.hh"
#include "paje.hh"
#include <getopt.h>
#include <string.h>

extern "C"
{
        extern FILE *yyin, *yyout;
        extern int yyline;
        int yyparse(void);
        int yylex(void);
        int yywrap()
        {
                return 1;
        }
}

using namespace std;


// ***********************************************************************
// Args parsing
//
struct _global_opts {
  string fout_name;
  string fin_name;
} global_opts;

static const char *optString = "i:o:h";

static const struct option longOpts[] = {
    { "description"           , required_argument , NULL, 'i' },
    { "output"                , required_argument , NULL, 'o' },
    { "help"                  , no_argument       , NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};

void usage(char ** argv)
{
  cout << "Usage: " <<argv[0] 
          << "<options>  file1.rst file2,rst ..."
          <<endl
       << "-i, --description  <file> : " 
           << "Description file for the rastro (default: stdin" 
           <<endl
       << "-o, --output <file>       : "
           << "Paje output from the rastro files (default: stdout)" 
           <<endl;
}

int parse_opts(int argc, char ** argv)
{
  int opt;
  int longIndex;

  global_opts.fout_name = "stdout";
  global_opts.fin_name  = "stdin";

  while (( opt = getopt_long( argc, argv, optString, longOpts, &longIndex ) ) != -1 )
  {
    switch( opt ) {
      case 'i':
          global_opts.fin_name = optarg;
          break;
      case 'o':
          global_opts.fout_name = optarg;
          break;
      case 'h':
          usage(argv);
          exit(0);
          break;
      case 0:
        if (!strcmp(longOpts[longIndex].name,"description"))
          global_opts.fin_name = optarg;
        else if (!strcmp(longOpts[longIndex].name,"output"))
          global_opts.fout_name = optarg;
        else if (!strcmp(longOpts[longIndex].name,"help"))
        {
          usage(argv);
          exit(0);
        }
        break;
      default:
        usage(argv);
        exit(1);
        break;
    }
  }

  return optind;
}



int main(int argc, char** argv)
{

  int oind;
  ostream *fout;
  FILE *fin; // dammed yacc

  oind = parse_opts(argc, argv);

  if (oind == argc) {
    cerr << "Expected some rst files as parameter" << endl;
    usage(argv);
    exit (1);
  }


  while (oind < argc) {
    cerr << "opening rastro " << argv[oind] <<endl;
    oind ++;
  }

  if (global_opts.fin_name == "stdin") {
    yyin = stdin;
  } else {
    fin = fopen(global_opts.fin_name.c_str(),"r");
    if (!fin) {
      cerr << "Could not open file " << global_opts.fin_name << endl;
      exit(1);
    }

    yyin = fin;
  }

  if (global_opts.fout_name == "stdout") {
    fout = &cout;
  } else {
    fout = new ofstream("out.test");
  }







  cout << "Hello World!" << endl;

  init_desc_parser();

  yyline = 1;
  yyparse();


  cout << "Container hierarchy at the end:" <<endl;
  print_tree(toplevel_hierarchy);


  check_unique_types();

  std::for_each(container_ids->begin(), container_ids->end(),
      [&](pair<string,hierarchy_t * > p) {
        hierarchy_t * cn = p.second;
        cout << p.first << " => " << *(cn->getVal()) ;
        if (cn->getParent())
          cout << " child of " << *(cn->getParent()->getVal()) ;
        cout << endl;
      } );


  init_paje_events();
  paje_header(*fout);
  hierarchy_to_paje(*fout);
  event_types_to_paje(*fout);

  ((ofstream *)fout)->close();

}
