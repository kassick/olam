// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Dom, 13 Nov 2011 02:11:40 -0200 (kassick)"
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
#include "paje.hh"
#include "rastro_loop.hh"
#include <getopt.h>
#include <string.h>
#include <set>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "rastro_generate.h"

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
    EXE_GEN_RST_FUNCTIONS_SOURCE
                          = 31,
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
  bool print_tree;


  // These are for the embeded rastro_generate
  bool has_rst_cheader, has_rst_csource, has_rst_fmodule;
  string rst_cheader, rst_csource, rst_fmodule;
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
    { "rst-c-functions"       , required_argument, NULL, 0 },
    { "rst-c-header"          , required_argument, NULL, 0 },
    { "rst-fort-module"       , required_argument, NULL, 0 },
    { "print-tree"            , no_argument      , NULL, 0 },
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
           <<endl
       << "--gen-ids <file>          : "
           << "Generate ID tokens for all events in the input file that do not have ids yet"
           <<endl
       << "--id-base <number>        : "
           << "Base for new IDs" 
           << endl
       << "--c-header <file>         : " 
           << "Generate a header file with defines for the event ids"
           <<endl
       << "--fort-header <file>      : "
           << "Generate a fortran module for the event ids"
           << endl
       << "--rst-signatures <file>   : "
           << "Generate a file compatible with rastro-generate"
           << endl
       << "--rst-c-functions         : "
           << "Generate a C file with the functions to call libRastro events"
           << endl
       << "--rst-c-header            : "
           << "Generate a C header file with the functions to call libRastro events"
           << endl
       << "--rst-fort-module         : "
           << "Generate a Fortran90 module file with the functions to call libRastro events"
           << endl
       << "--print-tree              : "
           << "Prints the tree after parsing the file"
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
  global_opts.print_tree = false;

  global_opts.has_rst_csource = false;
  global_opts.has_rst_cheader = false;
  global_opts.has_rst_fmodule = false;

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
        if (!strcmp(longOpts[longIndex].name,"description")) {
          global_opts.fin_name = optarg;
        } else if (!strcmp(longOpts[longIndex].name,"output")) {
          global_opts.fout_name = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "c-header")) {
          global_opts.exe_type.push_back(EXE_GEN_C_HEADER);
          global_opts.c_header = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "fort-header")) {
          global_opts.exe_type.push_back(EXE_GEN_FORT_HEADER);
          global_opts.fort_header = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "rst-signatures")) {
          global_opts.exe_type.push_back(EXE_GEN_RST_FUNCTION_SIGNATURES);
          global_opts.rst_signatures = optarg;
        } else if (!strcmp(longOpts[longIndex].name, "print-tree")) {
          global_opts.print_tree = true;
        } else if (!strcmp(longOpts[longIndex].name, "rst-c-functions")) {
          global_opts.has_rst_csource = true;
          global_opts.rst_csource = optarg;
          global_opts.exe_type.push_back(EXE_GEN_RST_FUNCTIONS_SOURCE);
        } else if (!strcmp(longOpts[longIndex].name, "rst-c-header")) {
          global_opts.has_rst_cheader = true;
          global_opts.rst_cheader = optarg;
          global_opts.exe_type.push_back(EXE_GEN_RST_FUNCTIONS_SOURCE);
        } else if (!strcmp(longOpts[longIndex].name, "rst-fort-module")) {
          global_opts.has_rst_fmodule = true;
          global_opts.rst_fmodule = optarg;
          global_opts.exe_type.push_back(EXE_GEN_RST_FUNCTIONS_SOURCE);
        } else if (!strcmp(longOpts[longIndex].name,"help")) {
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
      [&](pair<string, Paje::BaseEvent*> p) {
        Paje::BaseEvent *evt = p.second;

        if (!evt->has_ids()) {
          evt->gen_auto_ids(&(global_opts.auto_id_base), unique_ids);
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
  string module_name = filepath.stem().string();
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
      [&](pair<string, Paje::BaseEvent*> p) {
        Paje::BaseEvent *evt = p.second;
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
  string module_name = filepath.stem().string();

  h_file.open(fname);
  if (!h_file.good())
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }

  h_file << "Module " << module_name  << endl;


  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::BaseEvent*> p) {
        Paje::BaseEvent *evt = p.second;
          evt->gen_fort_header(h_file,"_IN","_OUT","_N");
      });
  
  
  h_file << "End Module " << module_name << endl;

  h_file.close();
}



/*
 *--------------------------------------------------------------------------------------
 *      Method: generate_rst_signatures_list
 *      Description:  Creates a list of signatures based on the events described
 *      in the input
 *--------------------------------------------------------------------------------------
 */
void generate_rst_signatures_list(set<string> & signatures) {
  for_each(ordered_event_names->begin(), ordered_event_names->end(),
      [&](pair<string, Paje::BaseEvent*> p) {
        p.second->get_rst_function_signature(signatures);
      });
}



void generate_rst_signatures_to_file(set<string> & signatures, FILE* fout) {
  for_each(signatures.begin(), signatures.end(),
      [&](string s) {
      if (s.size())
        fprintf(fout,"%s\n",s.c_str());
      } );
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  generate_rst_signatures
 *  Description:  Creates a file containing signatures to stand-alone
 *  rastro-generate
 * =====================================================================================
 */
void generate_rst_signatures(const string & fname)
{
  FILE *h_file;
  set<string> signatures;

  h_file  = fopen(fname.c_str(), "w");
  if (!h_file)
  {
    cerr << "[Error] Could not open file " << fname << endl;
    exit(1);
  }
  
  generate_rst_signatures_list(signatures);

  generate_rst_signatures_to_file(signatures,h_file);

  fclose(h_file);
}


FILE* generate_rst_signatures_to_tmpfile()
{
  FILE *h_file;
  set<string> signatures;

  h_file  = tmpfile();
  if (!h_file)
  {
    cerr << "[Error] Could not tmp file" << endl;
    exit(1);
  }
  
  generate_rst_signatures_list(signatures);

  generate_rst_signatures_to_file(signatures,h_file);
  rewind(h_file);

  return h_file;
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
  double timestamp = rastro_loop_events(rst_files_to_open,*fout,global_opts.print_tree); // loop every rst file and map events



  check_missing_links();

  //check_opened_states();

  check_opened_states(*fout);
  destroy_missing_containers(timestamp,*fout);
  
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
    files_to_parse.push(pathname.filename().string());
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


  if (global_opts.print_tree) {
    /// DEBUG!!!!
    print_tree(early_parse_tree);
    cerr << endl;
    print_tree(late_parse_tree);
    cerr << endl;
  }


  { // Do all the parsing of the files
    // Functions are in semantics.cc


    // Converts the early tree to the internal structures
    attr_to_container_hierarchy(early_parse_tree,toplevel_hierarchy);

    // create the types
    attr_to_event_types(early_parse_tree);
    attr_to_link_types(early_parse_tree);
    attr_to_state_types(early_parse_tree);

    // Create the events, links and dummy references
    attr_to_events(early_parse_tree);
    attr_to_states(early_parse_tree);
    attr_to_links(early_parse_tree);

    // parse the VALUE, EVENTTYPE and ID tokens that happen AFTER all the
    // definitions
    parse_late_tree();

    // Map types -> events/states/links
    map_accept_attrs(early_parse_tree);

    // add all events with ids to the ids map
    events_to_id_map();


    if (global_opts.print_tree) {
      // DEBUG DEBUG FEBUG DEBUH BUG BUG DE BUG

      cerr << "====================" << endl;
      cerr << "== Containers ======" << endl;
      for_each(container_type_names->begin(), container_type_names->end(),[&](pair<string,hierarchy_t * > p) {
          cerr << "container " << p.first << " defined" << endl;
        });

      cerr << "====================" << endl;
      cerr << "== Event Types =====" << endl;
      for_each(eventtype_names->begin(), eventtype_names->end(),[&](pair<string,Paje::BaseEventType * > p) {
          cerr << "event type " << p.first << " defined";
          if (p.second->container) 
            cerr << " on " << p.second->container->typeName;
          cerr << endl;
        });
      
      cerr << "====================" << endl;
      cerr << "== Event Names =====" << endl;
      for_each(event_names->begin(), event_names->end(), [&](pair<string,Paje::BaseEvent *> p) {
          cerr << p.second->toString() << endl;
          });

      cerr << "====================" << endl;
      cerr << "== Hierarchy =======" << endl;
      print_tree(toplevel_hierarchy);

      cerr << "====================" << endl;
      cerr << " == ID-Map =========" << endl;

      for_each(event_ids->begin(), event_ids->end(), [&]( pair<Paje::event_id_t, Paje::BaseEvent *> p) {
          cerr << p.first << " == > " << p.second->name << endl;
        });

    }
  } // Parsing




  sort(global_opts.exe_type.begin(),
       global_opts.exe_type.end() ); // luckly, there can be a hierarchy of what must be done before

  bool sources_generated = false;
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
          case EXE_GEN_RST_FUNCTIONS_SOURCE:
            if (!sources_generated)
            {
              FILE *tmp_signatures = generate_rst_signatures_to_tmpfile();
              generate_rst_functions(
                           global_opts.has_rst_csource, global_opts.rst_csource.c_str(),
                           global_opts.has_rst_cheader, global_opts.rst_cheader.c_str(),
                           global_opts.has_rst_fmodule, global_opts.rst_fmodule.c_str(),
                            tmp_signatures);
              sources_generated = true;
              fclose(tmp_signatures);
            }
            break;



          default:
            cerr << "Unknown exe type !? " << endl;
            exit(1);
        }
      });




  return 0;
}
