/* C source code
 * File: "/home/kassick/Work/olam/olamfs-trace-converter/src/olam_trace2paje.c"
 * Created: "Ter, 31 Mai 2011 11:11:38 -0300 (kassick)"
 * Updated: "Ter, 19 Jul 2011 19:26:43 -0300 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  olam_trace2paje.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31-05-2011 11:11:38 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include <search.h>
#include <rastro.h>
#include <assert.h>
#include <inttypes.h>

#include <aky.h>
#include <pvfs_events.h>
#include <olam_events.h>
#include <mpi_events.h>

#define FALSE (0)
#define TRUE (!FALSE)

#define OPTSTR "mho:v"
#define PROGNAME "olam_trace2paje"

#define OLAM_CONTAINER "APP"
#define PVFS_CONTAINER "FSPROCESS"
#define PAJE_ROOT_CONTAINER "0"

#define PROCESS_TYPE  "PROCESS"
#define FILE_TYPE     "FILE"
#define SERVER_TYPE   "SERVER"


void usage() {
  printf("Usage: %s <options> <file1.rst> <file2.rst> ...\n",PROGNAME);
  printf(" -h                 displays this message\n");
  printf(" -m                 adds MPI messages to the paje file\n");
  printf(" -v                 Adds the string values of some events to the value of paje events\n");
  printf(" -o <outputfile>    writes Pajé output to outputfile\n");
  //printf(" -o                 disables offset in constants (for olam rastro v>=11)");
  printf(" -somethingels      does some thing else\n");
}

int main(int argc, char** argv)
{
  int opt,nfiles = 0;
  char mpi_messages = FALSE;
  //int const_offset = OLAM_EVT_BASE;

  evt_name_t * ename;
  
  rst_file_t data;
  rst_event_t event;
  int i;
  int add_value_name = 0;

  //A function would do better...
  paje_ofile = stdout;

  while ((opt = getopt(argc,argv, OPTSTR)) != -1 )
  {
    switch(opt) {
      case 'm':
        mpi_messages = TRUE;
        break;

      case 'h':
        usage();
        exit(0);
        break;

      case 'o':
        paje_open_file(optarg);

        break;

      case 'v':
        add_value_name = TRUE;
        break;

        /*
      default:
        fprintf(stderr,"Unknown option -%c\n",(char)opt);
        usage();
        exit(1); */
    }
  }

  if ((argc - optind) == 0)
  {
    usage();
    exit(1);
  }


  // ---
  // Done with all the bureocracies


  hcreate(1000000); // Hash table to deal with lost arrows

  for (i = optind; i < argc; i++) {
    int ret = rst_open_file(argv[i], &data, NULL, 100000);
    if (ret == -1) {
      fprintf(stderr,"%s: trace %s could not be opened\n", argv[0], argv[i]);
      return 1;
    }
    nfiles++;
  }

  name_init();
  paje_header();
  paje_hierarchy();
  
  // Separate OLAM events (App and MPI Events) from FS events
  /*
  pajeCreateContainer(0,PVFS_CONTAINER,
                        PVFS_CONTAINER,
                        PAJE_ROOT_CONTAINER,
                        PVFS_CONTAINER);
*/
  while (rst_decode_event(&data, &event)) {
    char entity_name[100], app_name[100];
    char value[100];
    char state[100];
    char info[200];
    double timestamp;
    int type;
    char key[AKY_DEFAULT_STR_SIZE];
    char *locfn,*locfn_access;
    ENTRY he;
    

    if (event.id2 >= PVFS_VERSION_BASE)
    {
      // this is a PVFS event
      // Identified by server id and by file name
      snprintf(entity_name, 100, "pvfs_%" PRIu64 ,event.id1);
    } else { 
      // this is olam event -- identified by the process number
      snprintf(entity_name, 100, "rank%" PRIu64, event.id1); // this one is for container type PROCESS
      snprintf(app_name,    100, "olam_r%" PRIu64, event.id1); // this one for container type APP
    }
    
    snprintf(value, 100, "%s", name_get(event.type));
    ename = data_get(event.type); // only one of these will be valid, but ok
    timestamp = (double) event.timestamp / 1000000;

    if (add_value_name && (event.ct.n_string > 0)) {
      
      snprintf(info,200,"%s(%s",ename->short_name,event.v_string[0]);

      for (i = 1; i < event.ct.n_string; i++) {
        strncat(info,"_",200);
        strncat(info,event.v_string[i],200);
      }
      
      strcat(info,")");
      
      for (i = 0; i < strlen(info); i++)
      {
        if(info[i] == ' ') info[i] = '_';
      }
      //fprintf(stderr,"Adding value to event %s %s\n",ename->short_name,info);
    } else {
      snprintf(info,200,"%s",ename->short_name);
    }


    switch (event.type) {
      // OLAM and PVFS events here
      case OLAM_INIT:
      case MPI_INIT: // event_s -- hostname
        assert(event.ct.n_string == 1);

        //Container for the machine
        //Does it exist?
        he.key = event.v_string[0];
        if (!hsearch(he,FIND)) {
          if (!hsearch(he,ENTER)) {
            fprintf(stderr,"Could not add host %s to hash!\n",event.v_string[0]);
            exit(1);
          }
          fprintf(stderr,"got entry for machine %s\n",event.v_string[0]);
          pajeCreateContainer(timestamp,
                                event.v_string[0],
                                "MACHINE",
                                "0", event.v_string[0] );
        }
 
        // Create a container to the app
        pajeCreateContainer(timestamp,
                              app_name, // OLAM container
                              "APP", // of type APP
                              event.v_string[0],app_name); //in the machine

        // Create container to the rank
        pajeCreateContainer(timestamp,
                            entity_name, // rank_something
                            "PROCESS",
                            app_name, entity_name); //child of olam_something
        break;

      case PVFS_INIT:
        /*
        pajeCreateContainer(timestamp, entity_name,
                            SERVER_TYPE, PVFS_CONTAINER, entity_name);
                            */
        break;

      
#if 0
      case OLAM_HDF5_OPEN_IN:
        // Creates the container /filename/  -- Ignore for now
        // Creates arrow from rank%d to filename -- Ignore for now
        // Creates an state "opened" -- Ignore for now
        // open is iiss, myrank, thread_id, locfn,access
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        locfn_access = event.v_string[1];
     
        // Hierarchy has changed, there'll be no more duplicated file
        // elements

        //Create the container for the /locfn/ file
        // Type here is FILE -- Hope paje does not go crazy... ;)
        pajeCreateContainer(timestamp, locfn,
                            FILE_TYPE, PVFS_CONTAINER, locfn);
       
        //Push two stacked states to file:
        //  OPEN -- during the open operation
        //  OPENED -- from OPEN_IN to CLOSE_OUT
        snprintf(state,100, "STATE_%s",ename->start_name);
        pajePushState(timestamp, locfn, state, value);
        pajePushState(timestamp, locfn, "STATE_opened", "opened");

        // key for arrow:
        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"open %s#%" PRIu64 "#%s",
                                      locfn_access, event.id1, locfn);

        //Inter-container arrow -- ok?
        pajeStartLink(timestamp, OLAM_CONTAINER, "LINK", entity_name, "PTP", key);
        pajeEndLink  (timestamp, PVFS_CONTAINER, "LINK", locfn      , "PTP", key);
        
        break;
#endif


      case OLAM_SHDF5_INFO_IN:
      case OLAM_HDF5_CLOSE_READ_IN:
      case OLAM_SHDF5_OREC_IN:
      case OLAM_SHDF5_IREC_IN:
      case OLAM_HDF5_CLOSE_WRITE_IN:
      case OLAM_SHDF5_OPEN_IN:
      case OLAM_HDF5_DATASET_CLOSE_IN:
      case OLAM_HDF5_DATASET_GETINFO_IN:
      case OLAM_HDF5_DATASET_OPEN_IN:
      case OLAM_HDF5_PREPARE_READ_IN:
      case OLAM_HDF5_PREPARE_WRITE_IN:
      case OLAM_HDF5_OPEN_IN:
      case OLAM_HDF5_CREATE_IN:
      case OLAM_HDF5_WRITE_IN:
      case OLAM_HDF5_READ_IN:
        // Creates arrow from rank%d to filename
        // Creates an state with event

        //fprintf(stderr,"n = %d\n",event.ct.n_string);
        assert(event.ct.n_string >= 1);
        locfn = event.v_string[0];
        //locfn_access = event.v_string[1];
       
       // pop state open
        snprintf(state,100, "%s","P_STATE");
        //snprintf(state,100, "%s",ename->short_name);
        //snprintf(info,200,"%s %s",ename->short_name,locfn);

        pajePushState(timestamp, entity_name, state, info);

        /*
        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"%s#%" PRIu64 "#%s",
            entity_name, event.id1, locfn);
        pajeStartLink(timestamp, "0", "LINK", entity_name, "PTP", key);
        pajeEndLink  (timestamp, "0", "LINK", locfn      , "PTP", key);
        */
        break;

        /*
      case OLAM_HDF5_WRITE_IN:
      case OLAM_HDF5_READ_IN:
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        locfn_access = event.v_string[1];
      
        snprintf(state,100, "STATE_%s",ename->start_name);
        pajePushState(timestamp, locfn, state, value);

        // key for arrow:
        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"open %s#%" PRIu64 "#%s",
                                      locfn_access, event.id1, locfn);

        if (event.type == OLAM_HDF5_WRITE_IN) {
          pajeStartLink(timestamp, OLAM_CONTAINER, "LINK", entity_name, "PTP", key);
        } else {
          pajeStartLink(timestamp, PVFS_CONTAINER, "LINK", locfn,       "PTP", key);
        }


        break;

        */

        /*
      case OLAM_HDF5_WRITE_OUT:
      case OLAM_HDF5_READ_OUT:
      case OLAM_HDF5_OPEN_OUT:
      case OLAM_HDF5_CREATE_OUT:
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        locfn_access = event.v_string[1];
      
        snprintf(state,100, "STATE_%s",ename->start_name);
        pajePopState(timestamp, entity_name, state);

        // key for arrow:
        // id is "open W ## myrank ## /file/name"
        
        snprintf(key,AKY_DEFAULT_STR_SIZE,"open %s#%" PRIu64 "#%s",
                                      locfn_access, 
                                      event.id1, locfn);
        if (event.type == OLAM_HDF5_WRITE_OUT) {
          pajeEndLink  (timestamp, PVFS_CONTAINER, "LINK", locfn      , "PTP", key);
        } else {
          pajeEndLink  (timestamp, OLAM_CONTAINER, "LINK", entity_name, "PTP", key);
        }*/

        break;

      /*
      case OLAM_HDF5_CLOSE_OUT:
        // Destroys the container /filename/ 
        // Creates arrow from rank%d to filename
        // closes an state "opened"
        // open is iiss, myrank, thread_id, locfn,access
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        
        // pop close state and opened state
        snprintf(state,100, "STATE_%s",ename->start_name);
        pajePopState(timestamp, locfn, state);
        pajePopState(timestamp, locfn, "STATE_opened");
      
        pajeDestroyContainer(timestamp, FILE_TYPE, locfn);
        break;
      */


      case MPI_COMM_SPAWN_IN:
      case MPI_COMM_GET_NAME_IN:
      case MPI_COMM_SET_NAME_IN:
      case MPI_REDUCE_IN:
      case MPI_ALLREDUCE_IN:
      case MPI_REDUCE_SCATTER_IN:
      case MPI_ALLGATHER_IN:
      case MPI_ALLGATHERV_IN:
      case MPI_SCATTER_IN:
      case MPI_SCATTERV_IN:
      case MPI_WAIT_IN:
      case MPI_IRECV_IN:
      case MPI_ISEND_IN:
      case MPI_RECV_IN:
      case MPI_SEND_IN:
      case MPI_BCAST_IN:
      case MPI_BARRIER_IN:
      case MPI_GATHER_IN:
      case MPI_GATHERV_IN:
      case MPI_ALLTOALL_IN:
      case MPI_ALLTOALLV_IN:
      case MPI_OP_CREATE_IN:
      case MPI_OP_FREE_IN:
      case MPI_SCAN_IN:
      case MPI_ATTR_DELETE_IN:
      case MPI_ATTR_GET_IN:
      case MPI_ATTR_PUT_IN:
      case MPI_COMM_COMPARE_IN:
      case MPI_COMM_CREATE_IN:
      case MPI_COMM_DUP_IN:
      case MPI_COMM_FREE_IN:
      case MPI_COMM_GROUP_IN:
      case MPI_COMM_RANK_IN:
      case MPI_COMM_REMOTE_GROUP_IN:
      case MPI_COMM_REMOTE_SIZE_IN:
      case MPI_COMM_SIZE_IN:
      case MPI_COMM_SPLIT_IN:
      case MPI_COMM_TEST_INTER_IN:
      case MPI_GROUP_COMPARE_IN:
      case MPI_GROUP_DIFFERENCE_IN:
      case MPI_GROUP_EXCL_IN:
      case MPI_GROUP_FREE_IN:
      case MPI_GROUP_INCL_IN:
      case MPI_GROUP_INTERSECTION_IN:
      case MPI_GROUP_RANK_IN:
      case MPI_GROUP_RANGE_EXCL_IN:
      case MPI_GROUP_RANGE_INCL_IN:
      case MPI_GROUP_SIZE_IN:
      case MPI_GROUP_TRANSLATE_RANKS_IN:
      case MPI_GROUP_UNION_IN:
      case MPI_INTERCOMM_CREATE_IN:
      case MPI_INTERCOMM_MERGE_IN:
      case MPI_KEYVAL_CREATE_IN:
      case MPI_KEYVAL_FREE_IN:
      case MPI_ABORT_IN:
      case MPI_ERROR_CLASS_IN:
      case MPI_ERRHANDLER_CREATE_IN:
      case MPI_ERRHANDLER_FREE_IN:
      case MPI_ERRHANDLER_GET_IN:
      case MPI_ERROR_STRING_IN:
      case MPI_ERRHANDLER_SET_IN:
      case MPI_GET_PROCESSOR_NAME_IN:
      case MPI_INITIALIZED_IN:
      case MPI_WTICK_IN:
      case MPI_WTIME_IN:
      case MPI_ADDRESS_IN:
      case MPI_BSEND_IN:
      case MPI_BSEND_INIT_IN:
      case MPI_BUFFER_ATTACH_IN:
      case MPI_BUFFER_DETACH_IN:
      case MPI_CANCEL_IN:
      case MPI_REQUEST_FREE_IN:
      case MPI_RECV_INIT_IN:
      case MPI_SEND_INIT_IN:
      case MPI_GET_ELEMENTS_IN:
      case MPI_GET_COUNT_IN:
      case MPI_IBSEND_IN:
      case MPI_IPROBE_IN:
      case MPI_IRSEND_IN:
      case MPI_ISSEND_IN:
      case MPI_PACK_IN:
      case MPI_PACK_SIZE_IN:
      case MPI_PROBE_IN:
      case MPI_RSEND_IN:
      case MPI_RSEND_INIT_IN:
      case MPI_SENDRECV_IN:
      case MPI_SENDRECV_REPLACE_IN:
      case MPI_SSEND_IN:
      case MPI_SSEND_INIT_IN:
      case MPI_START_IN:
      case MPI_STARTALL_IN:
      case MPI_TEST_IN:
      case MPI_TESTALL_IN:
      case MPI_TESTANY_IN:
      case MPI_TEST_CANCELLED_IN:
      case MPI_TESTSOME_IN:
      case MPI_TYPE_COMMIT_IN:
      case MPI_TYPE_CONTIGUOUS_IN:
      case MPI_TYPE_EXTENT_IN:
      case MPI_TYPE_FREE_IN:
      case MPI_TYPE_HINDEXED_IN:
      case MPI_TYPE_HVECTOR_IN:
      case MPI_TYPE_INDEXED_IN:
      case MPI_TYPE_LB_IN:
      case MPI_TYPE_SIZE_IN:
      case MPI_TYPE_STRUCT_IN:
      case MPI_TYPE_UB_IN:
      case MPI_TYPE_VECTOR_IN:
      case MPI_UNPACK_IN:
      case MPI_WAITALL_IN:
      case MPI_WAITANY_IN:
      case MPI_WAITSOME_IN:
      case MPI_CART_COORDS_IN:
      case MPI_CART_CREATE_IN:
      case MPI_CART_GET_IN:
      case MPI_CART_MAP_IN:
      case MPI_CART_SHIFT_IN:
      case MPI_CARTDIM_GET_IN:
      case MPI_DIMS_CREATE_IN:
      case MPI_GRAPH_CREATE_IN:
      case MPI_GRAPH_GET_IN:
      case MPI_GRAPH_MAP_IN:
      case MPI_GRAPH_NEIGHBORS_IN:
      case MPI_GRAPH_NEIGHBORS_COUNT_IN:
      case MPI_GRAPHDIMS_GET_IN:
      case MPI_TOPO_TEST_IN:
      case MPI_RECV_IDLE_IN:
      case MPI_CART_RANK_IN:
      case MPI_CART_SUB_IN:
      case MPI_FINALIZE_IN:
        pajePushState(timestamp, entity_name, "MPI_STATE", value);
        break;

        /*
      case OLAM_HDF5_OPEN_IN: // create the state for open -- see later
      case OLAM_HDF5_CLOSE_IN:
      case OLAM_HDF5_CREATE_IN:
      case OLAM_HDF5_READ_IN:
      case OLAM_HDF5_WRITE_IN:
      */
      case OLAM_THREAD_IN:   // individual thread in/out
      case OLAM_PARBLOCK_IN: // whole parallel block
      case OLAM_READ_NL_IN:
      case OLAM_ONAME_CHECK_IN:
      case OLAM_HISTORY_START_IN:
      case OLAM_COPY_NL_IN:
      case OLAM_GRIDINIT_IN:
      case OLAM_PARA_DECOMP_IN:
      case OLAM_PARA_INIT_IN:
      case OLAM_MODSCHED_IN:
      case OLAM_FILL_JTABS_IN:
      case OLAM_FILL_JSEA_IN:
      case OLAM_FILL_JLAND_IN:
      case OLAM_FILL_JFLUX_IN:
      case OLAM_ALLOC_MISC_IN:
      case OLAM_JNMBINIT_IN:
      case OLAM_O_MEM_ALLOC_IN:
      case OLAM_O_ALLOC_MPI_IN:
      case OLAM_O_ALLOC_MPI_LAND_IN:
      case OLAM_O_ALLOC_MPI_SEA_IN:
      case OLAM_INITHH_IN:
      case OLAM_FLDSLHI_IN:  //WHAT!
      case OLAM_MICINIT_IN:
      case OLAM_HARR_RADINIT_IN:
      case OLAM_FULIOU_RADINIT_IN:
      case OLAM_LEAF3_STARTUP_IN:
      case OLAM_INIT_OFFLINE_MET_IN:
      case OLAM_READ_OFFLINE_MET_INIT_IN:
      case OLAM_LEAF3_INIT_ATM_IN:
      case OLAM_SEA_INIT_ATM_IN:
      case OLAM_ISAN_DRIVER_IN:
      case OLAM_RAYF_INIT_IN:
      case OLAM_PLOT_FIELDS_IN:
      case OLAM_HISTORY_WRITE_IN:
      case OLAM_MODEL_IN:
      case OLAM_O_CLSGKS_IN:
      case OLAM_TIMESTEP_IN:
      case OLAM_UPDATE_MODEL_TIME_IN:
      case OLAM_O_OUTPUT_IN:
      case OLAM_TEND0_IN:
      case OLAM_RADIATE_IN:
      case OLAM_SURFACE_TURB_FLUX_IN:
      case OLAM_TURB_K_IN:
      case OLAM_CUPARM_DRIVER_IN:
      case OLAM_SURFACE_CUPARM_FLUX_IN:
      case OLAM_THILTEND_LONG_IN:
      case OLAM_VELTEND_LONG_IN:
      case OLAM_OBS_NUDGE_IN:
      case OLAM_ZERO_MASSFLUX_IN:
      case OLAM_PROG_WRTU_IN:
      case OLAM_TIMEAVG_MASSFLUX_IN:
      case OLAM_SCALAR_TRANSPORT_IN:
      case OLAM_PREDTR_IN:
      case OLAM_THERMO_IN:
      case OLAM_MICRO_IN:
      case OLAM_SURFACE_PRECIP_FLUX_IN:
      case OLAM_TRSETS_IN:
      case OLAM_LEAF3_IN:
      case OLAM_SEACELLS_IN:
      case OLAM_INNERSTEP_IN:
        snprintf(state,100, "%s","APP_STATE");
        //snprintf(state,100, "%s",ename->short_name);
        pajePushState(timestamp, app_name, state,ename->short_name);

        break;

      //These ones use entity name as container
      case OLAM_OPLOT_INIT_IN:
      case OLAM_SEA_STARTUP_IN:
      case OLAM_SST_DATABASE_READ_IN:
      case OLAM_SEAICE_DATABASE_READ_IN:
      case OLAM_NDVI_DATABASE_READ_IN:
      case OLAM_ED_VEGETATION_DYNAMICS_IN:
      case OLAM_HDF5_CLOSE_IN:
      case OLAM_SHDF5_CLOSE_IN:
        snprintf(state,100, "%s","P_STATE");
        //snprintf(state,100, "%s",ename->short_name);
        pajePushState(timestamp, entity_name, state,ename->short_name);

        break;

      case MPI_COMM_SPAWN_OUT:
      case MPI_COMM_GET_NAME_OUT:
      case MPI_COMM_SET_NAME_OUT:
      case MPI_REDUCE_OUT:
      case MPI_ALLREDUCE_OUT:
      case MPI_REDUCE_SCATTER_OUT:
      case MPI_ALLGATHER_OUT:
      case MPI_ALLGATHERV_OUT:
      case MPI_SCATTER_OUT:
      case MPI_SCATTERV_OUT:
      case MPI_WAIT_OUT:
      case MPI_IRECV_OUT:
      case MPI_ISEND_OUT:
      case MPI_RECV_OUT:
      case MPI_SEND_OUT:
      case MPI_BCAST_OUT:
      case MPI_BARRIER_OUT:
      case MPI_GATHER_OUT:
      case MPI_GATHERV_OUT:
      case MPI_ALLTOALL_OUT:
      case MPI_ALLTOALLV_OUT:
      case MPI_OP_CREATE_OUT:
      case MPI_OP_FREE_OUT:
      case MPI_SCAN_OUT:
      case MPI_ATTR_DELETE_OUT:
      case MPI_ATTR_GET_OUT:
      case MPI_ATTR_PUT_OUT:
      case MPI_COMM_COMPARE_OUT:
      case MPI_COMM_CREATE_OUT:
      case MPI_COMM_DUP_OUT:
      case MPI_COMM_FREE_OUT:
      case MPI_COMM_GROUP_OUT:
      case MPI_COMM_RANK_OUT:
      case MPI_COMM_REMOTE_GROUP_OUT:
      case MPI_COMM_REMOTE_SIZE_OUT:
      case MPI_COMM_SIZE_OUT:
      case MPI_COMM_SPLIT_OUT:
      case MPI_COMM_TEST_INTER_OUT:
      case MPI_GROUP_COMPARE_OUT:
      case MPI_GROUP_DIFFERENCE_OUT:
      case MPI_GROUP_EXCL_OUT:
      case MPI_GROUP_FREE_OUT:
      case MPI_GROUP_INCL_OUT:
      case MPI_GROUP_INTERSECTION_OUT:
      case MPI_GROUP_RANK_OUT:
      case MPI_GROUP_RANGE_EXCL_OUT:
      case MPI_GROUP_RANGE_INCL_OUT:
      case MPI_GROUP_SIZE_OUT:
      case MPI_GROUP_TRANSLATE_RANKS_OUT:
      case MPI_GROUP_UNION_OUT:
      case MPI_INTERCOMM_CREATE_OUT:
      case MPI_INTERCOMM_MERGE_OUT:
      case MPI_KEYVAL_CREATE_OUT:
      case MPI_KEYVAL_FREE_OUT:
      case MPI_ABORT_OUT:
      case MPI_ERROR_CLASS_OUT:
      case MPI_ERRHANDLER_CREATE_OUT:
      case MPI_ERRHANDLER_FREE_OUT:
      case MPI_ERRHANDLER_GET_OUT:
      case MPI_ERROR_STRING_OUT:
      case MPI_ERRHANDLER_SET_OUT:
      case MPI_GET_PROCESSOR_NAME_OUT:
      case MPI_INITIALIZED_OUT:
      case MPI_WTICK_OUT:
      case MPI_WTIME_OUT:
      case MPI_ADDRESS_OUT:
      case MPI_BSEND_OUT:
      case MPI_BSEND_INIT_OUT:
      case MPI_BUFFER_ATTACH_OUT:
      case MPI_BUFFER_DETACH_OUT:
      case MPI_CANCEL_OUT:
      case MPI_REQUEST_FREE_OUT:
      case MPI_RECV_INIT_OUT:
      case MPI_SEND_INIT_OUT:
      case MPI_GET_ELEMENTS_OUT:
      case MPI_GET_COUNT_OUT:
      case MPI_IBSEND_OUT:
      case MPI_IPROBE_OUT:
      case MPI_IRSEND_OUT:
      case MPI_ISSEND_OUT:
      case MPI_PACK_OUT:
      case MPI_PACK_SIZE_OUT:
      case MPI_PROBE_OUT:
      case MPI_RSEND_OUT:
      case MPI_RSEND_INIT_OUT:
      case MPI_SENDRECV_OUT:
      case MPI_SENDRECV_REPLACE_OUT:
      case MPI_SSEND_OUT:
      case MPI_SSEND_INIT_OUT:
      case MPI_START_OUT:
      case MPI_STARTALL_OUT:
      case MPI_TEST_OUT:
      case MPI_TESTALL_OUT:
      case MPI_TESTANY_OUT:
      case MPI_TEST_CANCELLED_OUT:
      case MPI_TESTSOME_OUT:
      case MPI_TYPE_COMMIT_OUT:
      case MPI_TYPE_CONTIGUOUS_OUT:
      case MPI_TYPE_EXTENT_OUT:
      case MPI_TYPE_FREE_OUT:
      case MPI_TYPE_HINDEXED_OUT:
      case MPI_TYPE_HVECTOR_OUT:
      case MPI_TYPE_INDEXED_OUT:
      case MPI_TYPE_LB_OUT:
      case MPI_TYPE_SIZE_OUT:
      case MPI_TYPE_STRUCT_OUT:
      case MPI_TYPE_UB_OUT:
      case MPI_TYPE_VECTOR_OUT:
      case MPI_UNPACK_OUT:
      case MPI_WAITALL_OUT:
      case MPI_WAITANY_OUT:
      case MPI_WAITSOME_OUT:
      case MPI_CART_COORDS_OUT:
      case MPI_CART_CREATE_OUT:
      case MPI_CART_GET_OUT:
      case MPI_CART_MAP_OUT:
      case MPI_CART_SHIFT_OUT:
      case MPI_CARTDIM_GET_OUT:
      case MPI_DIMS_CREATE_OUT:
      case MPI_GRAPH_CREATE_OUT:
      case MPI_GRAPH_GET_OUT:
      case MPI_GRAPH_MAP_OUT:
      case MPI_GRAPH_NEIGHBORS_OUT:
      case MPI_GRAPH_NEIGHBORS_COUNT_OUT:
      case MPI_GRAPHDIMS_GET_OUT:
      case MPI_TOPO_TEST_OUT:
      case MPI_RECV_IDLE_OUT:
      case MPI_CART_RANK_OUT:
      case MPI_CART_SUB_OUT:
        pajePopState(timestamp, entity_name, "MPI_STATE");
        break;

        /*
      case OLAM_HDF5_OPEN_OUT:
      case OLAM_HDF5_CLOSE_OUT:
      case OLAM_HDF5_CREATE_OUT:
      case OLAM_HDF5_READ_OUT:
      case OLAM_HDF5_WRITE_OUT:
      */
      case OLAM_THREAD_OUT:   // individual thread in/out
      case OLAM_PARBLOCK_OUT: // whole parallel block
      case OLAM_READ_NL_OUT:
      case OLAM_ONAME_CHECK_OUT:
      case OLAM_HISTORY_START_OUT:
      case OLAM_COPY_NL_OUT:
      case OLAM_GRIDINIT_OUT:
      case OLAM_PARA_DECOMP_OUT:
      case OLAM_PARA_INIT_OUT:
      case OLAM_MODSCHED_OUT:
      case OLAM_FILL_JTABS_OUT:
      case OLAM_FILL_JSEA_OUT:
      case OLAM_FILL_JLAND_OUT:
      case OLAM_FILL_JFLUX_OUT:
      case OLAM_ALLOC_MISC_OUT:
      case OLAM_JNMBINIT_OUT:
      case OLAM_O_MEM_ALLOC_OUT:
      case OLAM_O_ALLOC_MPI_OUT:
      case OLAM_O_ALLOC_MPI_LAND_OUT:
      case OLAM_O_ALLOC_MPI_SEA_OUT:
      case OLAM_INITHH_OUT:
      case OLAM_FLDSLHI_OUT:  //WHAT!
      case OLAM_MICINIT_OUT:
      case OLAM_HARR_RADINIT_OUT:
      case OLAM_FULIOU_RADINIT_OUT:
      case OLAM_LEAF3_STARTUP_OUT:
      case OLAM_INIT_OFFLINE_MET_OUT:
      case OLAM_READ_OFFLINE_MET_INIT_OUT:
      case OLAM_LEAF3_INIT_ATM_OUT:
      case OLAM_SEA_INIT_ATM_OUT:
      case OLAM_ISAN_DRIVER_OUT:
      case OLAM_RAYF_INIT_OUT:
      case OLAM_PLOT_FIELDS_OUT:
      case OLAM_HISTORY_WRITE_OUT:
      case OLAM_MODEL_OUT:
      case OLAM_O_CLSGKS_OUT:
      case OLAM_TIMESTEP_OUT:
      case OLAM_UPDATE_MODEL_TIME_OUT:
      case OLAM_O_OUTPUT_OUT:
      case OLAM_TEND0_OUT:
      case OLAM_RADIATE_OUT:
      case OLAM_SURFACE_TURB_FLUX_OUT:
      case OLAM_TURB_K_OUT:
      case OLAM_CUPARM_DRIVER_OUT:
      case OLAM_SURFACE_CUPARM_FLUX_OUT:
      case OLAM_THILTEND_LONG_OUT:
      case OLAM_VELTEND_LONG_OUT:
      case OLAM_OBS_NUDGE_OUT:
      case OLAM_ZERO_MASSFLUX_OUT:
      case OLAM_PROG_WRTU_OUT:
      case OLAM_TIMEAVG_MASSFLUX_OUT:
      case OLAM_SCALAR_TRANSPORT_OUT:
      case OLAM_PREDTR_OUT:
      case OLAM_THERMO_OUT:
      case OLAM_MICRO_OUT:
      case OLAM_SURFACE_PRECIP_FLUX_OUT:
      case OLAM_TRSETS_OUT:
      case OLAM_LEAF3_OUT:
      case OLAM_SEACELLS_OUT:
      case OLAM_INNERSTEP_OUT:
        snprintf(state,100, "%s","APP_STATE");
        //snprintf(state,100, "%s",ename->short_name);
        pajePopState(timestamp, app_name, state);
        break;


      case OLAM_OPLOT_INIT_OUT:
      case OLAM_SEA_STARTUP_OUT:
      case OLAM_SST_DATABASE_READ_OUT:
      case OLAM_SEAICE_DATABASE_READ_OUT:
      case OLAM_NDVI_DATABASE_READ_OUT:
      case OLAM_ED_VEGETATION_DYNAMICS_OUT:
      case OLAM_SHDF5_CLOSE_OUT:

      case OLAM_SHDF5_INFO_OUT:
      case OLAM_HDF5_CLOSE_READ_OUT:
      case OLAM_SHDF5_OREC_OUT:
      case OLAM_SHDF5_IREC_OUT:
      case OLAM_HDF5_CLOSE_WRITE_OUT:
      case OLAM_SHDF5_OPEN_OUT:
      case OLAM_HDF5_DATASET_CLOSE_OUT:
      case OLAM_HDF5_DATASET_GETINFO_OUT:
      case OLAM_HDF5_DATASET_OPEN_OUT:
      case OLAM_HDF5_PREPARE_READ_OUT:
      case OLAM_HDF5_PREPARE_WRITE_OUT:
      case OLAM_HDF5_OPEN_OUT:
      case OLAM_HDF5_CREATE_OUT:
      case OLAM_HDF5_WRITE_OUT:
      case OLAM_HDF5_READ_OUT:
      case OLAM_HDF5_CLOSE_OUT:
        snprintf(state,100, "%s","P_STATE");
        //snprintf(state,100, "%s",ename->short_name);
        pajePopState(timestamp, entity_name, state);
        break;


      case MPI_FINALIZE_OUT:
        pajePopState(timestamp, entity_name, "MPI_STATE");
        pajeDestroyContainer(timestamp, "PROCESS", entity_name);
        pajeDestroyContainer(timestamp, "APP",     app_name);
        break;
      

#if 0
    case AKY_PTP_SEND:
      {
        char key[AKY_DEFAULT_STR_SIZE];
        aky_put_key("n", event.id1, event.v_uint32[0], key,
                    AKY_DEFAULT_STR_SIZE);
        pajeStartLink(timestamp, "0", "LINK", entity_name, "PTP", key);
      }
      break;
    case AKY_PTP_RECV:
      {
        char key[AKY_DEFAULT_STR_SIZE];
        aky_get_key("n", event.v_uint32[0], event.id1, key,
                    AKY_DEFAULT_STR_SIZE);
        pajeEndLink(timestamp, "0", "LINK", entity_name, "PTP", key);
      }
      break;
#endif
    }
  }

  rst_close_file(&data);
  hdestroy();
  return 0;

}
