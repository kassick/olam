// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/symbols.hh"
// Created: "Qui, 25 Ago 2011 14:38:26 -0300 (kassick)"
// Updated: "Seg, 29 Ago 2011 11:30:58 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  symbols.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  25-08-2011 14:38:26 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef __SYMBOLS_HH__
#define __SYMBOLS_HH__

#include <string>
#include <map>
#include <stdint.h>



using namespace std;

#define MAX_FORMAT_LEN 1024

#define decl_set_type(type,id) \
  void set_value(type arg)

#define impl_set_type(type,id) \
  void Symbol::set_value(type arg) { \
    this->val.id = arg;   \
    this->holds = id;  \
  }  \
  // trick?

namespace Paje {

  class Symbol {
    private:
      union {
        char c;
        uint16_t w;
        uint32_t i;
        uint64_t l;
        float f;
        double d;
        char *s;
      } val;

      enum {
        c, w, i, l, f, d, s,
      } holds;

    public:
      decl_set_type(char,c);
      decl_set_type(uint16_t,w);
      decl_set_type(uint32_t,i);
      decl_set_type(uint64_t,l);
      //decl_set_type(float,f);
      decl_set_type(double,d);
      decl_set_type(char*,s);

      void format(string fmt,stringstream &s);

  } ;

}




typedef map<string,Paje::Symbol *> symbols_table_t;


extern symbols_table_t * symbol_table;
void init_symbols();

string format_values(string & tpl, symbols_table_t *symbols);

#endif
