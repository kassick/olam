// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/symbols.hh"
// Created: "Qui, 25 Ago 2011 14:38:26 -0300 (kassick)"
// Updated: "Qui, 20 Out 2011 22:42:23 -0200 (kassick)"
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

#include <iostream>
#include <string>
#include <map>
#include <stdint.h>

#include "rastro_helper.hh"



using namespace std;

#define MAX_FORMAT_LEN 1024

#define decl_set_type(type,id) \
  void set_value(type arg)


#define SYMBOL_FREE_STR_PRE \
    if (this->holds == s) { \
      free (this->val.s);   \
      this->holds = invalid;\
    }

#define SYMBOL_NO_PRE
#define SYMBOL_NO_FUNCT

#define impl_set_type_funct(type,id,pre,funct) \
  void Symbol::set_value(type arg) { \
    pre \
    this->val.id = funct(arg);   \
    this->holds = id;  \
  }  \
  // trick?

#define impl_set_type(type,id) impl_set_type_funct(type, id, SYMBOL_FREE_STR_PRE, SYMBOL_NO_FUNCT)

namespace Paje {

  extern string idf1_name, idf2_name;
  

  class Symbol {
    private:
      rastro_basic_val_t val;

      rastro_basic_types_t holds;


    public:
      decl_set_type(uint8_t,c);
      decl_set_type(char,c);
      decl_set_type(uint16_t,w);
      decl_set_type(uint32_t,i);
      decl_set_type(uint64_t,l);
      //decl_set_type(float,f);
      decl_set_type(double,d);
      decl_set_type(const char*,s);

      Symbol& operator= (const Symbol & org);
      void format(const string &fmt,ostream &s);

      Symbol();
      Symbol(const Symbol &s);
      ~Symbol();

  } ;

}



extern bool symbols_format_ok;


typedef map<string,Paje::Symbol> symbols_table_t;


extern symbols_table_t * symbol_table;
void init_symbols();

string format_values(const string & tpl, symbols_table_t **symbols,bool warn = true);
bool   format_values(const string & tpl, symbols_table_t **symbols,ostream &out,
                            bool warn = true);

#endif
