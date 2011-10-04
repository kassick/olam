// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/paje_functions.cc"
// Created: "Ter, 04 Out 2011 14:14:35 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:17:48 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  paje_functions.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 14:14:35 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include "paje_functions.hh"

using namespace std;

// Helper functions adapted from akypuera
static double first_timestamp = -1;

//FILE * paje_ofile = NULL;


// These are private to this file; no need to be on the header
typedef struct paje_event {
  string name;
  string description;
  int id;

  paje_event(const string& _name, const string& _description) {
    name = _name;
    description = _description;
    id = -1;
  }
} paje_event_t;


typedef map<string,paje_event_t *> paje_events_t;

paje_events_t * paje_events;



static int paje_event_id(const char *name)
{
  paje_events_t::iterator it;
  it = paje_events->find(name);
  if (it == paje_events->end())
      return -1;

  return (*it).second->id;
}

static double paje_event_timestamp(double timestamp)
{
  if (first_timestamp == -1) {
    first_timestamp = timestamp;
  }
  return timestamp - first_timestamp;
}


//============================================================
// Public Functions:
//
void init_paje_events() {
  paje_events = new paje_events_t();
  (*paje_events)["PajeDefineContainerType"] = 
    new paje_event_t("PajeDefineContainerType",
       "% Alias string\n"
       "% ContainerType string\n"
       "% Name string");
  (*paje_events)["PajeDefineStateType"] = 
    new paje_event_t("PajeDefineStateType",
       "% Alias string\n"
       "% ContainerType string\n"
       "% Name string");
  (*paje_events)["PajeDefineLinkType"] = 
    new paje_event_t("PajeDefineLinkType",
       "% Alias string\n"
       "% ContainerType string\n"
       "% SourceContainerType string\n"
       "% DestContainerType string\n"
       "% Name string");
  (*paje_events)["PajeCreateContainer"] = 
    new paje_event_t("PajeCreateContainer",
       "% Time string\n"
       "% Alias string\n"
       "% Type string\n"
       "% Container string\n"
       "% Name string");
  (*paje_events)["PajeDestroyContainer"] = 
    new paje_event_t("PajeDestroyContainer",
       "% Time string\n"
       "% Type string\n"
       "% Container string");
  (*paje_events)["PajeSetState"] =
    new paje_event_t("PajeSetState",
       "% Time string\n"
       "% Container string\n"
       "% EntityType string\n"
       "% Value string");
  (*paje_events)["PajePushState"] = 
    new paje_event_t("PajePushState",
       "% Time string\n"
       "% Container string\n"
       "% EntityType string\n"
       "% Value string");
  (*paje_events)["PajePopState"] = 
    new paje_event_t("PajePopState",
       "% Time string\n"
       "% Container string\n"
       "% EntityType string");
  (*paje_events)["PajeStartLink"] = 
    new paje_event_t("PajeStartLink",
       "% Time string\n"
       "% Container string\n"
       "% EntityType string\n"
       "% SourceContainer string\n"
       "% Value string\n"
       "% Key string");
  (*paje_events)["PajeEndLink"] = 
    new paje_event_t("PajeEndLink",
       "% Time string\n"
       "% Container string\n"
       "% EntityType string\n"
       "% DestContainer string\n"
       "% Value string\n"
       "% Key string" );

  //number the events
  paje_events_t::iterator it;
  int id = 0;
  for (it = paje_events->begin(); it != paje_events->end(); ++it)
  {
    paje_event_t * evt = (*it).second;
    evt->id = id++;
  }
};

void pajeDefineContainerType(string &alias,
                             string &containerType, 
                             string &name,
                             ostream &out)
{
  out << paje_event_id("PajeDefineContainerType")
      << " " << alias
      << " " << containerType
      << " " << name
      <<endl;
}

void pajeDefineStateType(string &alias,
                         string &containerType, 
                         string &name,
                         ostream &out)
      
{
  out << paje_event_id("PajeDefineStateType")
      << " " << alias
      << " " << containerType
      << " " << name
      << endl;
}

void pajeDefineLinkType(string &alias,
                        string &containerType,
                        string &sourceContainerType,
                        string &destContainerType, 
                        string &name,
                        ostream &out)
{
  out << paje_event_id("PajeDefineLinkType")
      << " " << alias
      << " " << containerType
      << " " << sourceContainerType
      << " " << destContainerType
      << " " <<   name
      <<endl;
}

void pajeCreateContainer(double timestamp,
                         string &alias,
                         string &type,
                         string &container, string &name,
                         ostream &out)
{
  out << paje_event_id("PajeCreateContainer")
      << " " << paje_event_timestamp(timestamp)
      << " " << alias
      << " " << type
      << " " << container
      << " " << name
      << endl;
}

void pajeDestroyContainer(double timestamp,
                          string &type, string &container,
                          ostream &out)
{
  out << paje_event_id("PajeDestroyContainer")
      << " " << paje_event_timestamp(timestamp)
      << " " << type
      << " " << container
      << endl;
}

void pajeSetState(double timestamp,
                  string &container,
                  string &type, string &value,
                  ostream &out)
{
  out  << paje_event_id("PajeSetState")
       << " " << paje_event_timestamp(timestamp)
       << " " << container
       << " " << type
       << " " << value
       << endl;
}

void pajePushState(double timestamp,
                   string &container,
                   string &type, string &value,
                   ostream &out)
{
  out << paje_event_id("PajePushState")
      << " " << paje_event_timestamp(timestamp)
      << " " << container
      << " " << type
      << " " << value
      << endl;
}

void pajePopState(double timestamp,
                  string &container, string &type,
                  ostream &out)
{
  out << paje_event_id("PajePopState")
      << " " << paje_event_timestamp(timestamp)
      << " " << container
      << " " << type
      << endl;
}

void pajeStartLink(double timestamp,
                   string &container,
                   string &type,
                   string &sourceContainer,
                   string &value, string &key,
                   ostream &out)
{
  out << paje_event_id("PajeStartLink")
      << " " << paje_event_timestamp(timestamp)
      << " " << container
      << " " << type
      << " " << sourceContainer
      << " " << value
      << " " << key
      << endl;
}

void pajeEndLink(double timestamp,
                 string &container,
                 string &type,
                 string &endContainer,
                 string &value, string &key,
                 ostream &out)
{
  out << paje_event_id("PajeEndLink")
      << " " << paje_event_timestamp(timestamp)
      << " " << container
      << " " << type
      << " " << endContainer
      << " " << value
      << " " << key
      << endl;
}

void paje_header(ostream &out)
{
  paje_events_t::iterator it;
  for (it = paje_events->begin(); it != paje_events->end(); ++it)
  {
    //cout << "evt " <<(*it).first << " is " << (*it).second <<endl;
    paje_event_t * evt = (*it).second;
    out << "%EventDef " << evt->name << " " << evt->id <<endl;
    out << evt->description << endl;
    out << "%EndEventDef" <<endl;
  }
}

#if 0
#define STATE_NAME_MAX 300
void paje_hierarchy(void)
{
  int i;
  pajeDefineContainerType("MACHINE",    "0"         , "MACHINE");
  pajeDefineContainerType("APP",        "MACHINE"   , "APP");
  pajeDefineContainerType("FILESYSTEM", "MACHINE"   , "FS");
  pajeDefineContainerType("OMP_THREAD", "APP"       , "OMP");
  pajeDefineContainerType("PROCESS"   , "APP"       , "PROC");
  pajeDefineContainerType("FILE"      , "FILESYSTEM", "FILE");
  pajeDefineContainerType("FSPROCESS" , "FILESYSTEM", "FSPROC");

  pajeDefineStateType("APP_STATE"    , "APP"    , "APP_STATE"); // App has app wide states
  pajeDefineStateType("P_STATE"    , "PROCESS"    , "P_STATE"); // App has app wide states
  pajeDefineStateType("T_STATE"    , "OMP_THREAD"    , "T_STATE"); // App has app wide states
  
  //pajeDefineStateType("APP_STATE"    , "APP"    , "APP_STATE"); // App has app wide states
  pajeDefineStateType("MPI_STATE", "PROCESS", "MPI_STATE"); // App has app wide states
  // Now add an state type for each event in olam
  /*
  for (i = 0; olam_evt_names[i].name != NULL; i++) {
    char state_name[STATE_NAME_MAX];
    snprintf(state_name,STATE_NAME_MAX,"%s",olam_evt_names[i].short_name);
        //start_name);
    pajeDefineStateType(state_name, "PROCESS", state_name);
    pajeDefineStateType(state_name, "APP", state_name);
  }
  */

  /*
  for (i = 0; pvfs_evt_names[i].name != NULL; i++) {
    char state_name[STATE_NAME_MAX];
    snprintf(state_name,STATE_NAME_MAX,"STATE_%s",pvfs_evt_names[i].start_name);
    pajeDefineStateType(state_name, "FSPROCESS", state_name);
    pajeDefineStateType(state_name, "FILE"     , state_name);
  }
  */
  
  pajeDefineStateType("STATE_OPENED", "FILE"   , "STATE_OPENED");

  //pajeDefineLinkType("COMMLINK"  , "APP"    , "PROCESS", "PROCESS" , "COMMLINK");
  //pajeDefineLinkType("FILELINK"  , "MACHINE", "PROCESS", "FILE"    , "FILELINK");
  //pajeDefineLinkType("FILEOPLINK", "MACHINE", "PROCESS", "FSPROC"  , "FILEOPLINK");
}



int paje_open_file(const char const * fname)
{
  paje_ofile = fopen(fname,"w");
  if (!paje_ofile)
  {
    fprintf(stderr,"Cannot open %s for writing, fallback to stdout\n",fname);
    paje_ofile = stdout;
    return 1;
  }
  return 0;
}

#endif
