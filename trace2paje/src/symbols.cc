// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/symbols.cc"
// Created: "Qui, 25 Ago 2011 14:03:15 -0300 (kassick)"
// Updated: "Seg, 29 Ago 2011 12:13:16 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  symbols.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  25-08-2011 14:03:15 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include "symbols.hh"
#include <stdio.h>
#include <string.h>
#include <cinttypes>
#include <boost/xpressive/xpressive.hpp>


using namespace std;
using namespace boost::xpressive;

symbols_table_t * symbol_table;

void init_symbols() {
  symbol_table = new symbols_table_t();
}


namespace Paje {

  impl_set_type(char,c)
  impl_set_type(uint16_t,w)
  impl_set_type(uint32_t,i)
  impl_set_type(uint64_t,l)
  //impl_set_type(float,f)
  impl_set_type(double,d)
  impl_set_type(char*,s)

  void Symbol::format(string fmt,stringstream &ss)
  {
    string format =  "%" + fmt;
    char cc = fmt[fmt.length()-1];
    char tmpbuf[MAX_FORMAT_LEN];

    if (fmt.length() == 0 || !('a' <= cc <= 'z')) 
    {
      // user has not specified a format; it's up to us!
      switch(this->holds) 
      {
        case f:
          format += "f";
          break;
        case d:
          format += "f";
          break;
        case c:
          format += "c";
          break;
        case w:
          format += PRIu16;
          break;
        case i:
          format += PRIu32;
          break;
        case l:
          format += PRIu64;
          break;
        case s:
          format += "s";
          break;
      }

    }

    snprintf(tmpbuf,MAX_FORMAT_LEN,format.c_str(),this->val);

    ss << tmpbuf;

  }


}


#define ID_REGEX "[a-zA-Z][a-zA-Z0-9_]*"
#define FMT_REGEX "[#0-9\\.\\- \\+'IhlLqjztdiouxXeEfFgGaAcsCSpnm]+"
//#define FMT_REGEX "[#0-9]+"

string format_values(string & tpl, symbols_table_t *symbols)
{
  stringstream out;
  string tmp = tpl;
  smatch res;
  string fmt_regex("%\\((?P<id>" ID_REGEX ")(:(?P<fmt>" FMT_REGEX "))?\\)");
  //string fmt_regex("%\\((P<id>[a-zA-Z][a-zA-Z0-9_]*)(:(P?<fmt>.+))?\\)");
  sregex rx = sregex::compile(fmt_regex);

  while (tmp.length() > 0) {
    string id;
    string fmt;
    if (regex_search(tmp,res,rx)) {

      out << res.prefix().str();

      try{
        id = res["id"];
      } catch (...) {
        // should not happen, id is a non optional arg
        cerr <<"no id for you" <<endl;
        exit(1);
      }

      try {
        fmt = res["fmt"];
      } catch (...) {
        // no format, let default
        fmt = "";
      }


      symbols_table_t::iterator it;
      it = symbols->find(id);
      if (it != symbols->end()) {
        Paje::Symbol *s = (*it).second;
        s->format(fmt,out);
      } else {
        cerr << "No symbol ``" << id << "´´, ignoring" << endl;
      }
    }

    tmp = res.suffix().str();
  }

  return out.str();

}

int main(int argc, char ** argv)
{
  cout << "hello" <<endl;
  Paje::Symbol * s1;
  stringstream s;
  string id;

  init_symbols();
  
  s1 = new Paje::Symbol();
  s1->set_value((uint32_t)10);
  id = "int1";
  (*symbol_table)[id] = s1;
  
  s1 = new Paje::Symbol();
  s1->set_value("hello");
  (*symbol_table)[string("str1")] = s1;
  
  s1 = new Paje::Symbol();
  s1->set_value((float)56.23453);
  (*symbol_table)[string("float1")] = s1;
  
  s1 = new Paje::Symbol();
  s1->set_value((double)52.23);
  (*symbol_table)[string("double1")] = s1;


  string tpl = "this is a prefix of %(null:10) %(str1:.10s) is followed by %(float1) and %(double1:10.2G)";

  cout << format_values(tpl,symbol_table) << endl;
 
  tpl = "nothing to see here";
  cout << format_values(tpl,symbol_table) << endl;



  /*
  s << "inicio ";
  

  s1->format("20.3d",s);
  s << " " ;
  s1->format("20.3o",s);

  s1->format("",s);
  s1->format("20s",s);
  s << endl;

  s1->format("",s);
  s << endl;

  cout << s.str();

  cout << (float)56.23453 <<endl;
  */

  return 0;
}

