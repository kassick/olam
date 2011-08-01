// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Seg, 01 Ago 2011 18:02:40 -0300 (kassick)"
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

#include <iostream>
#include <algorithm>
#include "semantics.hh"
#include "container.hh"

extern "C"
{
        extern int yyline;
        int yyparse(void);
        int yylex(void);
        int yywrap()
        {
                return 1;
        }
}

using namespace std;

int main(int argc, char** argv)
{
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


}
