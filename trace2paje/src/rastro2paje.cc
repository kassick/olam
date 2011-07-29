// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro2paje.cc"
// Created: "Ter, 26 Jul 2011 13:01:06 -0300 (kassick)"
// Updated: "Qui, 28 Jul 2011 19:53:50 -0300 (kassick)"
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
#include "container.hh"

extern "C"
{
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

  /* This is test code */
  attribs_t *test_tree;
  struct semantic_attribute * attr;
  attribs_t::iterator it1,it2;
  test_tree = new attribs_t();
  attr = new semantic_attribute();
  attr->id = ID_NOP;

  it1 = test_tree->insert(test_tree->begin(),attr);
  
  attr = new semantic_attribute();
  attr->id = ID_NOP;
  it2 = test_tree->insert(test_tree->begin(),attr);

  cout <<"flat tree" <<endl;
  print_tree(test_tree,test_tree->begin(),test_tree->end());

  test_tree->reparent(it2,it1);

  cout <<"just it2" <<endl;
  print_tree(test_tree,it2,it2.end());

  cout << "whole tree" << endl;
  print_tree(test_tree,test_tree->begin(),test_tree->end());

  cout <<"out!" <<endl;

  exit(0);

  


  init_desc_parser();

  yyparse();


}
