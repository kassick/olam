// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/symbols.cc"
// Created: "Qui, 25 Ago 2011 14:03:15 -0300 (kassick)"
// Updated: "Qui, 25 Ago 2011 18:11:04 -0300 (kassick)"
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

using namespace std;

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



int main(int argc, char ** argv)
{
  cout << "hello" <<endl;
  Paje::Symbol * s1;
  stringstream s;


  s << "inicio ";
  
  s1 = new Paje::Symbol();
  s1->set_value((uint32_t)10);

  s1->format("20.3d",s);
  s << " " ;
  s1->format("20.3o",s);

  s1->set_value("hello");
  s1->format("",s);
  s1->format("20s",s);
  s << endl;

  s1->set_value((float)56.23453);
  s1->format("",s);
  s << endl;

  cout << s.str();

  cout << (float)56.23453 <<endl;

  return 0;
}

