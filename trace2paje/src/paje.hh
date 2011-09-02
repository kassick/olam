// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/paje.hh"
// Created: "Seg, 01 Ago 2011 15:34:40 -0300 (kassick)"
// Updated: "Sex, 02 Set 2011 17:51:58 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  paje.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  01-08-2011 15:34:40 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef _PAJE_HH_H
#define _PAJE_HH_H

#include <string>
#include <sstream>

using namespace std;

namespace Paje {

  class PajeElement {
    public:
      void toPaje(stringstream &definitions, stringstream &header, stringstream &finalization);
  };
}

void init_paje_events();
void paje_header(ostream &out);
void pajeDefineContainerType(string &alias,
                             string &containerType, 
                             string &name,
                             ostream &out);

void pajeDefineStateType(string &alias,
                         string &containerType, 
                         string &name,
                         ostream &out);

void pajeDefineLinkType(string &alias,
                        string &containerType,
                        string &sourceContainerType,
                        string &destContainerType, string &name,
                        ostream &out);

void pajeCreateContainer(double timestamp,
                         string &alias,
                         string &type,
                         string &container, string &name,
                         ostream &out);

void pajeDestroyContainer(double timestamp,
                          string &type, string &container,
                          ostream &out);

void pajeSetState(double timestamp,
                  string &container,
                  string &type, string &value,
                  ostream &out);

void pajePushState(double timestamp,
                   string &container,
                   string &type, string &value,
                   ostream &out);

void pajePopState(double timestamp,
                  string &container, string &type,
                  ostream &out);

void pajeStartLink(double timestamp,
                   string &container,
                   string &type,
                   string &sourceContainer,
                   string &value, string &key,
                   ostream &out);

void pajeEndLink(double timestamp,
                 string &container,
                 string &type,
                 string &endContainer,
                 string &value, string &key,
                 ostream &out);

#endif
