/*
    This file is part of Akypuera

    Akypuera is free software: you can redistribute it and/or modify
    it under the terms of the GNU Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Akypuera is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Public License for more details.

    You should have received a copy of the GNU Public License
    along with Akypuera. If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <aky.h>
#include <rastro.h>

//static timestamp_t first_timestamp = -1;
static double first_timestamp = -1;
FILE * paje_ofile = NULL;

static s_paje_event_t paje_events[] = {
  {"PajeDefineContainerType",
   "% Alias string\n" "% ContainerType string\n" "% Name string",
   -1},
  {"PajeDefineStateType",
   "% Alias string\n" "% ContainerType string\n" "% Name string",
   -1},
  {"PajeDefineLinkType",
   "% Alias string\n"
   "% ContainerType string\n"
   "% SourceContainerType string\n"
   "% DestContainerType string\n" "% Name string",
   -1},
  {"PajeCreateContainer",
   "% Time string\n"
   "% Alias string\n"
   "% Type string\n" "% Container string\n" "% Name string",
   -1},
  {"PajeDestroyContainer",
   "% Time string\n" "% Type string\n" "% Container string",
   -1},
  {"PajeSetState",
   "% Time string\n"
   "% Container string\n" "% EntityType string\n" "% Value string",
   -1},
  {"PajePushState",
   "% Time string\n"
   "% Container string\n" "% EntityType string\n" "% Value string",
   -1},
  {"PajePopState",
   "% Time string\n" "% Container string\n" "% EntityType string",
   -1},
  {"PajeStartLink",
   "% Time string\n"
   "% Container string\n"
   "% EntityType string\n"
   "% SourceContainer string\n" "% Value string\n" "% Key string",
   -1},
  {"PajeEndLink",
   "% Time string\n"
   "% Container string\n"
   "% EntityType string\n"
   "% DestContainer string\n" "% Value string\n" "% Key string",
   -1},
  {NULL, NULL, -1}
};

static int paje_event_id(const char *name)
{
  int i;
  for (i = 0; paje_events[i].name; i++)
    if (strcmp(name, paje_events[i].name) == 0)
      return paje_events[i].id;
  return -1;
}

static double paje_event_timestamp(double timestamp)
{
  if (first_timestamp == -1) {
    first_timestamp = timestamp;
  }
  return timestamp - first_timestamp;
}

void pajeDefineContainerType(const char *alias,
                             const char *containerType, const char *name)
{
  fprintf(paje_ofile,"%d %s %s %s\n",
         paje_event_id("PajeDefineContainerType"),
         alias, containerType, name);
}

void pajeDefineStateType(const char *alias,
                         const char *containerType, const char *name)
{
  fprintf(paje_ofile,"%d %s %s %s\n",
         paje_event_id("PajeDefineStateType"), alias, containerType, name);
}

void pajeDefineLinkType(const char *alias,
                        const char *containerType,
                        const char *sourceContainerType,
                        const char *destContainerType, const char *name)
{
  fprintf(paje_ofile,"%d %s %s %s %s %s\n",
         paje_event_id("PajeDefineLinkType"),
         alias, containerType, sourceContainerType, destContainerType,
         name);
}

void pajeCreateContainer(double timestamp,
                         const char *alias,
                         const char *type,
                         const char *container, const char *name)
{
  fprintf(paje_ofile,"%d %f %s %s %s %s\n",
         paje_event_id("PajeCreateContainer"),
         paje_event_timestamp(timestamp), alias, type, container, name);
}

void pajeDestroyContainer(double timestamp,
                          const char *type, const char *container)
{
  fprintf(paje_ofile,"%d %f %s %s\n",
         paje_event_id("PajeDestroyContainer"),
         paje_event_timestamp(timestamp), type, container);
}

void pajeSetState(double timestamp,
                  const char *container,
                  const char *type, const char *value)
{
  fprintf(paje_ofile,"%d %f %s %s %s\n",
         paje_event_id("PajeSetState"),
         paje_event_timestamp(timestamp), container, type, value);
}

void pajePushState(double timestamp,
                   const char *container,
                   const char *type, const char *value)
{
  fprintf(paje_ofile,"%d %f %s %s %s\n",
         paje_event_id("PajePushState"),
         paje_event_timestamp(timestamp), container, type, value);
}

void pajePopState(double timestamp,
                  const char *container, const char *type)
{
  fprintf(paje_ofile,"%d %f %s %s\n",
         paje_event_id("PajePopState"),
         paje_event_timestamp(timestamp), container, type);
}

void pajeStartLink(double timestamp,
                   const char *container,
                   const char *type,
                   const char *sourceContainer,
                   const char *value, const char *key)
{
  fprintf(paje_ofile,"%d %f %s %s %s %s %s\n",
         paje_event_id("PajeStartLink"),
         paje_event_timestamp(timestamp),
         container, type, sourceContainer, value, key);
}

void pajeEndLink(double timestamp,
                 const char *container,
                 const char *type,
                 const char *endContainer,
                 const char *value, const char *key)
{
  fprintf(paje_ofile,"%d %f %s %s %s %s %s\n",
         paje_event_id("PajeEndLink"),
         paje_event_timestamp(timestamp),
         container, type, endContainer, value, key);
}

void paje_header(void)
{
  int i;
  for (i = 0; paje_events[i].name; i++) {
    paje_events[i].id = i;
    fprintf(paje_ofile,"%%EventDef %s %d\n%s\n%%EndEventDef\n",
           paje_events[i].name, paje_events[i].id,
           paje_events[i].description);
  }
}

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
