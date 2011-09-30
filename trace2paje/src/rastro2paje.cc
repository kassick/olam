// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Sex, 30 Set 2011 19:30:48 -0300 (kassick)"
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
#include <queue>
#include "semantics.hh"
#include "container.hh"
#include "paje.hh"
#include "rastro_loop.hh"
#include <getopt.h>
#include <string.h>
#include <set>
#include <boost/filesystem.hpp>

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
namespace fs = boost::filesystem;


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
  
  fs::path pathname; // path of the first file
  fs::path dirname; // dirname of the first file

  list<string> rst_files_to_open;

  oind = parse_opts(argc, argv);

  if (oind == argc) {
    cerr << "Expected some rst files as parameter" << endl;
    usage(argv);
    exit (1);
  }


  while (oind < argc) {
    rst_files_to_open.push_back(string(argv[oind]));
    oind ++;
  }


  if (global_opts.fin_name == "stdin") {
    pathname = "";
    files_to_parse.push("/dev/stdin");
    //yyin = stdin;
  } else {
    pathname = global_opts.fin_name;
    files_to_parse.push(pathname.filename());
    /*
    fin = fopen(global_opts.fin_name.c_str(),"r");
    if (!fin) {
      cerr << "Could not open file " << global_opts.fin_name << endl;
      exit(1);
    }

    yyin = fin;
    */
  }

  init_desc_parser();
  // Here we parse every description file we have in the INCLUDE directives

  set<string> parsed_files;
 

  while (!files_to_parse.empty()) {
    string file_name = files_to_parse.front();
    files_to_parse.pop();

    dirname = pathname.parent_path();
    dirname /= file_name;

    string fin_name = dirname.string();


    if (parsed_files.find(fin_name) != parsed_files.end())
    {
      cerr << "Ignoring file " << fin_name << " (already parsed)" << endl;
    } else {

      cerr << "Parsing "<<fin_name << "... ";
      
      fin = fopen(fin_name.c_str(),"r");
      if (!fin) {
        cerr << "Could not open file " << global_opts.fin_name << endl;
        exit(1);
      }

      yyin = fin;
      yyline = 1;
      yyparse();

      cout << " [done]" << endl;
      fclose(fin);

      parsed_files.insert(fin_name);
    
    
    }



  }


#if 0
  /// DEBUG!!!!
  print_tree(early_parse_tree);
  cerr << endl;
  print_tree(late_parse_tree);
  cerr << endl;
  /*
  walk_tree_depth_first(late_parse_tree,[&](attribs_t * n, int level) {
      while (level--)
        cout << " ";
      cout << "attr id is " << n->getVal()->id;
      cout << endl;
      return false;
    });
    */
#endif


  // Converts the early tree to the internal structures
  attr_to_container_hierarchy(early_parse_tree,toplevel_hierarchy);

  // create the types
  attr_to_event_types(early_parse_tree);
  attr_to_link_types(early_parse_tree);



#if 0
  for_each(container_type_names->begin(), container_type_names->end(),[&](pair<string,hierarchy_t * > p) {
      cerr << "container " << p.first << " defined" << endl;
    });
#endif


#if 0
  for_each(eventtype_names->begin(), eventtype_names->end(),[&](pair<string,Paje::EventType * > p) {
      cerr << "event " << p.first << " defined" << endl;
    });

  cerr << "Container hierarchy at the end:" <<endl;
  print_tree(toplevel_hierarchy);

#endif

  // Create the events, links and dummy references

  attr_to_states(early_parse_tree);
  attr_to_links(early_parse_tree);


  // parse the VALUE, EVENTTYPE and ID tokens that happen AFTER all the
  // definitions
  parse_late_tree();

  // Map types -> events/states/links
  map_accept_attrs(early_parse_tree);

  // add all events with ids to the ids map
  events_to_id_map();


#if 0
  for_each(event_names->begin(), event_names->end(), [&](pair<string,Paje::Event *> p) {
      cerr << p.second->toString() << endl;
      });
#endif


  // Here begins the parsing of events and the output

  if (global_opts.fout_name == "stdout") {
    fout = &cout;
  } else {
    fout = new ofstream(global_opts.fout_name);
    if (fout->fail())
    {
      cerr << "[Error] Can not open file " << global_opts.fout_name << endl;
      exit(1);
    }
  }


  /*
  std::for_each(container_type_names->begin(), container_type_names->end(),
      [&](pair<string,hierarchy_t * > p) {
        hierarchy_t * cn = p.second;
        cerr << p.first << " => " << *(cn->getVal()) ;
        if (cn->getParent())
          cerr << " child of " << *(cn->getParent()->getVal()) ;
        cerr << endl;
      } ); */



  // Here's the output
  init_paje_events();
  paje_header(*fout);
  hierarchy_to_paje(*fout);
  event_types_to_paje(*fout);

  rastro_loop_events(rst_files_to_open,*fout);


  // Here loop around the paje events and dumps the events

  ((ofstream *)fout)->close();

}
