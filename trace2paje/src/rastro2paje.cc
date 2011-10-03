// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Seg, 03 Out 2011 18:14:29 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  rastro2paje.cc
 *
 *    Description:  Take a description file and convert the given rastros
 *                  to Pajé Language
 *
 *        Version:  1.0
 *        Created:  26-07-2011 13:01:06 BRT
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:   (r v kassick), Rodrigo Virote Kassick
 *        Company:  GPPD/II-UFRGS
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
#include <boost/algorithm/string.hpp>

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

typedef enum _EXE_TYPE {
    EXE_GEN_IDS           = 0,
    EXE_GEN_C_HEADER      = 10,
    EXE_GEN_FORT_HEADER   = 20,
    EXE_GEN_RST_FUNCTION_SIGNATURES 
                          = 30,
    EXE_GEN_PAJE          = 100,
    EXE_NOP               = 101,
  } exe_enum_t;

// ***********************************************************************
// Args parsing
//
struct _global_opts {
  string fout_name;
  string fin_name;
  vector<exe_enum_t> exe_type;
  long int auto_id_base;
  string ids_file, c_header, fort_header, rst_signatures;
} global_opts;

static const char *optString = "i:o:hg:";

static const struct option longOpts[] = {
    { "description"           , required_argument, NULL, 'i' },
    { "output"                , required_argument, NULL, 'o' },
    { "help"                  , no_argument      , NULL, 'h' },
    { "gen-ids"               , required_argument, NULL, 'g' },
    { "id-base"               , required_argument, NULL, 'b' },
    { "c-header"              , required_argument, NULL, 0 },
    { "fort-header"           , required_argument, NULL, 0 },
    { "rst-signatures"        , required_argument, NULL, 0 },
    { NULL, no_argument, NULL, 0  }
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

#define AUTO_ID_BASE_DEFAULT 4000

int parse_opts(int argc, char ** argv)
{
  int opt;
  int longIndex;

  global_opts.fout_name = "stdout";
  global_opts.fin_name  = "stdin";
  global_opts.auto_id_base = AUTO_ID_BASE_DEFAULT;

  while (( opt = getopt_long( argc, argv, optString, longOpts, &longIndex ) ) != -1 )
  {
    switch( opt ) {
      case 'i':
          global_opts.fin_name = optarg;
          break;
      case 'o':
          global_opts.fout_name = optarg;
          global_opts.exe_type.push_back(EXE_GEN_PAJE);
          break;
      case 'h':
          usage(argv);
          exit(0);
          break;
      case 'g':
          global_opts.exe_type.push_back(EXE_GEN_IDS);
          global_opts.ids_file = optarg;
          break;
      case 'b':
          global_opts.auto_id_base = atol(optarg);
          break;
      case 0:
        if (!strcmp(longOpts[longIndex].name,"description"))
          global_opts.fin_name = optarg;
        else if (!strcmp(longOpts[longIndex].name,"output"))
          global_opts.fout_name = optarg;
        else if (!strcmp(longOpts[longIndex].name, "c-header")) {
          global_opts.exe_type.push_back(EXE_GEN_C_HEADER);
          global_opts.c_header = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "fort-header")) {
          global_opts.exe_type.push_back(EXE_GEN_FORT_HEADER);
          global_opts.fort_header = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "rst-signatures")) {
          global_opts.exe_type.push_back(EXE_GEN_RST_FUNCTION_SIGNATURES);
          global_opts.rst_signatures = optarg;
        }
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

  if (global_opts.exe_type.size() == 0) {
      global_opts.exe_type.push_back(EXE_GEN_PAJE);
  }

  return optind;
}


void generate_ids_to_file(const string & fname)
{
  ofstream ids_file;

  ids_file.open(fname);

  if (!ids_file.good())
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }

  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::Event*> p) {
        Paje::Event *evt = p.second;

        if (!evt->has_ids()) {
          evt->gen_auto_ids(&(global_opts.auto_id_base));
          evt->gen_ids_description(ids_file);
        }

      });
  ids_file.close();
}



void generate_c_header(const string & fname)
{
  ofstream h_file;
  fs::path filepath;
  filepath = fname;
  string module_name = filepath.stem();
  boost::to_upper(module_name);

  h_file.open(fname);

  if (!h_file.good())
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }

  h_file << "#ifndef __" << module_name << "_H__" << endl
        <<  "#define __" << module_name << "_H__" << endl;
  /*  we do not need extern C if not defining functions
        <<  "#ifdef __cplusplus"    << endl
        <<  "extern \"C\" {"        << endl
        <<  "#endif"                << endl;*/


  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::Event*> p) {
        Paje::Event *evt = p.second;
          evt->gen_c_header(h_file,"_IN","_OUT","_N");
      });
  
  
  /*h_file << "#ifdef __cplusplus"      << endl
         << "}"                       << endl
         << "#endif"                  << endl*/
  h_file << "#endif  //__"<< module_name << "_H__" << endl;

  h_file.close();
}

void generate_fort_header(const string & fname)
{
  fs::path filepath;
  ofstream h_file;

  filepath = fname;
  string module_name = filepath.stem();

  h_file.open(fname);
  if (!h_file.good())
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }

  h_file << "Module " << module_name  << endl;


  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::Event*> p) {
        Paje::Event *evt = p.second;
          evt->gen_fort_header(h_file,"_IN","_OUT","_N");
      });
  
  
  h_file << "End Module " << module_name << endl;

  h_file.close();
}

void generate_rst_signatures(const string & fname)
{
  ofstream h_file;
  set<string> signatures;

  h_file.open(fname);
  if (!h_file.good())
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }
  
  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::Event*> p) {
        p.second->get_rst_function_signature(signatures);
      });

  for_each(signatures.begin(), signatures.end(),
      [&](string s) {
        h_file << s << endl;
      } );

  h_file.close();
}

void generate_paje_output(const string & fname,int oind, int argc, char ** argv)
{
  ostream *fout;
  list<string> rst_files_to_open;       // .rst files
  // Set output file
  
  if (fname == "stdout") {
    fout = &cout;
  } else {
    fout = new ofstream(fname);
    if (fout->fail())
    {
      cerr << "[Error] Can not open file " << global_opts.fout_name << endl;
      exit(1);
    }
  }

  // get the rst files from the argument line
  if (oind == argc) {
    cerr << "Expected some rst files as parameter" << endl;
    usage(argv);
    exit (1);
  }

  while (oind < argc) {
    rst_files_to_open.push_back(string(argv[oind]));
    oind ++;
  }

  // OUTPUT:
  init_paje_events();         // INIT paje functions (helper)
  paje_header(*fout);         // Define Paje functions in the file
  hierarchy_to_paje(*fout);   // Define Container Hierarchy
  event_types_to_paje(*fout); // Define event types

  // Here loop around the paje events and dumps the events
  rastro_loop_events(rst_files_to_open,*fout); // loop every rst file and map events

  ((ofstream *)fout)->close();

}


int main(int argc, char** argv)
{
  int oind;
  FILE *fin; // dammed yacc
  
  fs::path pathname; // path of the first file
  fs::path dirname; // dirname of the first file
  // files to parse is defined in semantic_types.hh / semantics.cc -- it's
  // needed globally to parse the INPUT statements


  oind = parse_opts(argc, argv);


  
  //default to stdin
  if (global_opts.fin_name == "stdin") {
    pathname = "";
    files_to_parse.push("/dev/stdin");
    //yyin = stdin;
  } else {
    pathname = global_opts.fin_name;
    files_to_parse.push(pathname.filename());
  }

  // Creates basic structures
  init_desc_parser();

  
  
  // Here we parse every description file we have in the INCLUDE directives
  set<string> parsed_files;
  while (!files_to_parse.empty()) {
    string file_name = files_to_parse.front();
    files_to_parse.pop();

    // all files are relative to the first path
    dirname = pathname.parent_path();
    dirname /= file_name;

    string fin_name = dirname.string();


    // parse each file only once
    if (parsed_files.count(fin_name))
    {
      cerr << "Ignoring file " << fin_name << " (already parsed)" << endl;
    } else {

      cerr << "Parsing "<<fin_name << "... ";
      
      // lex/yacc depend on standard C file handles :(
      fin = fopen(fin_name.c_str(),"r");
      if (!fin) {
        cerr << "Could not open file " << global_opts.fin_name << endl;
        exit(1);
      }

      yyin = fin;
      yyline = 1;
      yyparse(); // Yacc!!

      cout << " [done]" << endl;
      fclose(fin);

      parsed_files.insert(fin_name);
    }
  } // while (!files_to_parse.empty())


#if 0
  /// DEBUG!!!!
  print_tree(early_parse_tree);
  cerr << endl;
  print_tree(late_parse_tree);
  cerr << endl;
#endif


  { // Do all the parsing of the files
    // Functions are in semantics.cc


    // Converts the early tree to the internal structures
    attr_to_container_hierarchy(early_parse_tree,toplevel_hierarchy);

    // create the types
    attr_to_event_types(early_parse_tree);
    attr_to_link_types(early_parse_tree);

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
    // DEBUG DEBUG FEBUG DEBUH BUG BUG DE BUG
    for_each(container_type_names->begin(), container_type_names->end(),[&](pair<string,hierarchy_t * > p) {
        cerr << "container " << p.first << " defined" << endl;
      });

    for_each(eventtype_names->begin(), eventtype_names->end(),[&](pair<string,Paje::EventType * > p) {
        cerr << "event " << p.first << " defined" << endl;
      });

    cerr << "Container hierarchy at the end:" <<endl;
    print_tree(toplevel_hierarchy);

    for_each(event_names->begin(), event_names->end(), [&](pair<string,Paje::Event *> p) {
        cerr << p.second->toString() << endl;
        });
#endif
  } // Parsing




  sort(global_opts.exe_type.begin(),
       global_opts.exe_type.end() ); // luckly, there can be a hierarchy of what must be done before

  for_each(global_opts.exe_type.begin(),global_opts.exe_type.end(),
      [&](exe_enum_t exe_type) {

        switch (exe_type) {
          case EXE_GEN_IDS:
            generate_ids_to_file(global_opts.ids_file);
            break;
          case EXE_GEN_PAJE:
            generate_paje_output(global_opts.fout_name,oind, argc, argv);
            break;
          case EXE_GEN_C_HEADER:
            generate_c_header(global_opts.c_header);
            break;
          case EXE_GEN_FORT_HEADER:
            generate_fort_header(global_opts.fort_header);
            break;
          case EXE_GEN_RST_FUNCTION_SIGNATURES:
            generate_rst_signatures(global_opts.rst_signatures);
            break;


          default:
            cerr << "Unknown exe type !? " << endl;
            exit(1);
        }
      });




  return 0;
}
