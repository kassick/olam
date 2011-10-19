// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/symbols.cc"
// Created: "Qui, 25 Ago 2011 14:03:15 -0300 (kassick)"
// Updated: "Qua, 19 Out 2011 17:03:54 -0200 (kassick)"
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

  string idf1_name, idf2_name;

  impl_set_type(uint8_t,c)
  impl_set_type(char,c)
  impl_set_type(uint16_t,w)
  impl_set_type(uint32_t,i)
  impl_set_type(uint64_t,l)
  //impl_set_type(float,f)
  impl_set_type(double,d)
  //impl_set_type_funct(const char*,s, SYMBOL_FREE_STR_PRE, strdup)
  void Symbol::set_value(const char * str1) {
    if (this->holds == s)
      free(this->val.s);

    this->val.s = strdup(str1);
    this->holds = s;
  }


  Symbol::Symbol() {
    //cout <<" Hey, created a new one at" << this << endl;
    this->holds = invalid;
  }

  Symbol::Symbol(const Symbol& org) {
    this->holds = org.holds;
    if (this->holds == s)
      this->val.s == strdup(org.val.s);
    else
      memcpy( & (this->val), & (org.val), sizeof(org.val));
  }


  Symbol::~Symbol() {
    //SYMBOL_FREE_STR_PRE ;
    if (this->holds == s)
      free(this->val.s);
  }

  Symbol& Symbol::operator= (const Symbol & org)
  {
    if (this == &org)
      return *this;

    if (this->holds == s)
      free(this->val.s);

    this->holds = org.holds;
    if (this->holds == s)
      this->val.s = strdup(org.val.s);
    else
      memcpy( & (this->val), & (org.val), sizeof(org.val));

    return *this;
  }


  void Symbol::format(string fmt,ostream &ss)
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


bool symbols_format_ok;

string format_values(string & tpl,
    symbols_table_t **symbols,
    bool warn) {
  stringstream out;

  symbols_format_ok = format_values(tpl,symbols,out,warn);

  return out.str();
}


bool format_values(string & tpl, 
    symbols_table_t **symbols,
    ostream &out, bool warn)
{
  bool all_ok = true;
  string tmp = tpl;
  smatch res;
  string fmt_regex("%\\((?P<id>" ID_REGEX ")(:(?P<fmt>" FMT_REGEX "))?\\)");
  //string fmt_regex("%\\((P<id>[a-zA-Z][a-zA-Z0-9_]*)(:(P?<fmt>.+))?\\)");
  sregex rx = sregex::compile(fmt_regex);

  while (tmp.length() > 0) {
    string id;
    string fmt;
    bool had_fmt = false;

    bool found = regex_search(tmp,res,rx);
    /*
    cerr << "+ prefix: " << res.prefix().str() <<endl;
    cerr << "+ suffix: " << res.suffix().str() <<endl;
    cerr << "+ tmp:    " << tmp <<endl;
    */

    if (found) {

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
        had_fmt = true;
      } catch (...) {
        // no format, let default
        fmt = "";
      }


      symbols_table_t::iterator it;
      int i;
      for (i = 0;symbols[i] != NULL;i++)
      {
        if ((it = symbols[i]->find(id)) != symbols[i]->end())
        {
          Paje::Symbol &s = (*it).second;
          s.format(fmt,out);
          break;
        }
      }

      if (symbols[i] == NULL)
      {
        // unknown symbol; just dump the fmt string 
        out << "(" << id;
        if (had_fmt)
          out << ":" << fmt;
        out << ")";
        if (warn)
          cerr << "No symbol ``" << id << "´´, ignoring" << endl;

        all_ok = false;
      }
    
      tmp = res.suffix().str();
    } else {
      //cerr << "not found" <<endl;
      out << tmp;
      tmp = "";
    }

  }

  return all_ok;

}


#ifdef SYMBOL_UNIT_TEST
int main(int argc, char ** argv)
{
  cout << "hello" <<endl;
  Paje::Symbol * s1, s2, s3;
  stringstream s;
  string id;

  cout << "init symbols and create s1" << endl;
  init_symbols();
 
  s1 = new Paje::Symbol();



  cout << "Symbol uint32" <<endl;
  s1->set_value((uint32_t)10);
  id = "int1";
  (*symbol_table)[id] = *s1;
  
  cout << "Symbol str" <<endl;
  s1->set_value("hello");
  (*symbol_table)[string("str1")] = *s1;
  
  cout << "Symbol float" <<endl;
  s1->set_value((float)56.23453);
  (*symbol_table)[string("float1")] = *s1;
  
  cout << "Symbol double" <<endl;
  s1->set_value((double)52.23);
  (*symbol_table)[string("double1")] = *s1;


  cout << "now format it" << endl;


  string tpl = "this is a prefix of %(null:10) %(str1:.10s) is followed by %(float1) and %(double1:e)";

  symbols_table_t* v[2];
  v[0] = symbol_table;
  v[1] = NULL;

  cout << "template: " << tpl << endl;
  cout << format_values(tpl,v) << endl;
  cout << "---" << endl;
  
  tpl = "this is a prefix of %(null:10) %(str1:.10s) is followed by %(float1) and %(double1:10.2G) and a trailiing";
  cout << "template: " << tpl << endl;
  cout << format_values(tpl,v) << endl;
  cout << "---" << endl;
 
  tpl = "nothing to see here";
  cout << "template: " << tpl << endl;
  cout << format_values(tpl,v) << endl;
  cout << "---" << endl;



  cout << "Now test copy operator " << endl;

  s2.set_value("hello motto");

  s3 = s2;

  (*symbol_table)["s2"] = s2;
  (*symbol_table)["s3"] = s3;

  Paje::Symbol &tmp = (*symbol_table)["undefined"];
  tmp.set_value("undefined I am");

  tpl = "s2 is %(s2) and s3 is %(s3) \"%(undefined)\"";

  format_values(tpl, v, cout);
  cout << endl;



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

#endif
