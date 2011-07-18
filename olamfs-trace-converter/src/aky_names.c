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
#include <string.h>
#include <stdlib.h>
#include "aky.h"
#include <olam_events.h>
#include <pvfs_events.h>
#include <mpi_events.h>



static aky_name_entry_t * mpi_names = NULL;




void name_init(void)
{
  int i;
  mpi_names = malloc(sizeof(aky_name_entry_t) * MAX_AKY_STATE_NAMES);
  bzero(mpi_names,sizeof(aky_name_entry_t) * MAX_AKY_STATE_NAMES);
#if 0

  mpi_names[MPI_INIT].name = strdup("MPI_Init");
  mpi_names[MPI_COMM_SPAWN_IN].name = strdup("MPI_Comm_spawn");
  mpi_names[MPI_COMM_GET_NAME_IN].name = strdup("MPI_Comm_get_name");
  mpi_names[MPI_COMM_SET_NAME_IN].name = strdup("MPI_Comm_set_name");
  mpi_names[MPI_REDUCE_IN].name = strdup("MPI_Reduce");
  mpi_names[MPI_ALLREDUCE_IN].name = strdup("MPI_Allreduce");
  mpi_names[MPI_REDUCE_SCATTER_IN].name = strdup("MPI_Reduce_scatter");
  mpi_names[MPI_ALLGATHER_IN].name = strdup("MPI_Allgather");
  mpi_names[MPI_ALLGATHERV_IN].name = strdup("MPI_Allgatherv");
  mpi_names[MPI_SCATTER_IN].name = strdup("MPI_Scatter");
  mpi_names[MPI_SCATTERV_IN].name = strdup("MPI_Scatterv");
  mpi_names[MPI_WAIT_IN].name = strdup("MPI_Wait");
  mpi_names[MPI_IRECV_IN].name = strdup("MPI_Irecv");
  mpi_names[MPI_ISEND_IN].name = strdup("MPI_Isend");
  mpi_names[MPI_RECV_IN].name = strdup("MPI_Recv");
  mpi_names[MPI_SEND_IN].name = strdup("MPI_Send");
  mpi_names[MPI_BCAST_IN].name = strdup("MPI_Bcast");
  mpi_names[MPI_BARRIER_IN].name = strdup("MPI_Barrier");
  mpi_names[MPI_GATHER_IN].name = strdup("MPI_Gather");
  mpi_names[MPI_GATHERV_IN].name = strdup("MPI_Gatherv");
  mpi_names[MPI_ALLTOALL_IN].name = strdup("MPI_Alltoall");
  mpi_names[MPI_ALLTOALLV_IN].name = strdup("MPI_Alltoallv");
  mpi_names[MPI_OP_CREATE_IN].name = strdup("MPI_Op_create");
  mpi_names[MPI_OP_FREE_IN].name = strdup("MPI_Op_free");
  mpi_names[MPI_SCAN_IN].name = strdup("MPI_Scan");
  mpi_names[MPI_ATTR_DELETE_IN].name = strdup("MPI_Attr_delete");
  mpi_names[MPI_ATTR_GET_IN].name = strdup("MPI_Attr_get");
  mpi_names[MPI_ATTR_PUT_IN].name = strdup("MPI_Attr_put");
  mpi_names[MPI_COMM_COMPARE_IN].name = strdup("MPI_Comm_compare");
  mpi_names[MPI_COMM_CREATE_IN].name = strdup("MPI_Comm_create");
  mpi_names[MPI_COMM_DUP_IN].name = strdup("MPI_Comm_dup");
  mpi_names[MPI_COMM_FREE_IN].name = strdup("MPI_Comm_free");
  mpi_names[MPI_COMM_GROUP_IN].name = strdup("MPI_Comm_group");
  mpi_names[MPI_COMM_RANK_IN].name = strdup("MPI_Comm_rank");
  mpi_names[MPI_COMM_REMOTE_GROUP_IN].name = strdup("MPI_Comm_remote_group");
  mpi_names[MPI_COMM_REMOTE_SIZE_IN].name = strdup("MPI_Comm_remote_size");
  mpi_names[MPI_COMM_SIZE_IN].name = strdup("MPI_Comm_size");
  mpi_names[MPI_COMM_SPLIT_IN].name = strdup("MPI_Comm_split");
  mpi_names[MPI_COMM_TEST_INTER_IN].name = strdup("MPI_Comm_test_inter");
  mpi_names[MPI_COMM_TEST_INTER_OUT].name = strdup("MPI_Comm_test_inter_out");
  mpi_names[MPI_GROUP_COMPARE_IN].name = strdup("MPI_Group_compare");
  mpi_names[MPI_GROUP_DIFFERENCE_IN].name = strdup("MPI_Group_difference");
  mpi_names[MPI_GROUP_EXCL_IN].name = strdup("MPI_Group_excl");
  mpi_names[MPI_GROUP_FREE_IN].name = strdup("MPI_Group_free");
  mpi_names[MPI_GROUP_INCL_IN].name = strdup("MPI_Group_incl");
  mpi_names[MPI_GROUP_INCL_OUT].name = strdup("MPI_Group_incl_out");
  mpi_names[MPI_GROUP_INTERSECTION_IN].name = strdup("MPI_Group_intersection");
  mpi_names[MPI_GROUP_INTERSECTION_OUT].name = strdup("MPI_Group_intersection_out");
  mpi_names[MPI_GROUP_RANK_IN].name = strdup("MPI_Group_rank");
  mpi_names[MPI_GROUP_RANGE_EXCL_IN].name = strdup("MPI_Group_range_excl");
  mpi_names[MPI_GROUP_RANGE_INCL_IN].name = strdup("MPI_Group_range_incl");
  mpi_names[MPI_GROUP_RANGE_INCL_OUT].name = strdup("MPI_Group_range_incl_out");
  mpi_names[MPI_GROUP_SIZE_IN].name = strdup("MPI_Group_size");
  mpi_names[MPI_GROUP_TRANSLATE_RANKS_IN].name = strdup("MPI_Group_translate_ranks");
  mpi_names[MPI_GROUP_UNION_IN].name = strdup("MPI_Group_union");
  mpi_names[MPI_INTERCOMM_CREATE_IN].name = strdup("MPI_Intercomm_create");
  mpi_names[MPI_INTERCOMM_CREATE_OUT].name = strdup("MPI_Intercomm_create_out");
  mpi_names[MPI_INTERCOMM_MERGE_IN].name = strdup("MPI_Intercomm_merge");
  mpi_names[MPI_INTERCOMM_MERGE_OUT].name = strdup("MPI_Intercomm_merge_out");
  mpi_names[MPI_KEYVAL_CREATE_IN].name = strdup("MPI_Keyval_create");
  mpi_names[MPI_KEYVAL_FREE_IN].name = strdup("MPI_Keyval_free");
  mpi_names[MPI_ABORT_IN].name = strdup("MPI_Abort");
  mpi_names[MPI_ERROR_CLASS_IN].name = strdup("MPI_Error_class");
  mpi_names[MPI_ERRHANDLER_CREATE_IN].name = strdup("MPI_Errhandler_create");
  mpi_names[MPI_ERRHANDLER_FREE_IN].name = strdup("MPI_Errhandler_free");
  mpi_names[MPI_ERRHANDLER_GET_IN].name = strdup("MPI_Errhandler_get");
  mpi_names[MPI_ERROR_STRING_IN].name = strdup("MPI_Error_string");
  mpi_names[MPI_ERRHANDLER_SET_IN].name = strdup("MPI_Errhandler_set");
  mpi_names[MPI_GET_PROCESSOR_NAME_IN].name = strdup("MPI_Get_processor_name");
  mpi_names[MPI_INITIALIZED_IN].name = strdup("MPI_Initialized");
  mpi_names[MPI_INITIALIZED_OUT].name = strdup("MPI_Initialized_out");
  mpi_names[MPI_WTICK_IN].name = strdup("MPI_Wtick");
  mpi_names[MPI_WTIME_IN].name = strdup("MPI_Wtime");
  mpi_names[MPI_ADDRESS_IN].name = strdup("MPI_Address");
  mpi_names[MPI_BSEND_IN].name = strdup("MPI_Bsend");
  mpi_names[MPI_BSEND_INIT_IN].name = strdup("MPI_Bsend_init");
  mpi_names[MPI_BSEND_INIT_OUT].name = strdup("MPI_Bsend_init_out");
  mpi_names[MPI_BUFFER_ATTACH_IN].name = strdup("MPI_Buffer_attach");
  mpi_names[MPI_BUFFER_DETACH_IN].name = strdup("MPI_Buffer_detach");
  mpi_names[MPI_CANCEL_IN].name = strdup("MPI_Cancel");
  mpi_names[MPI_REQUEST_FREE_IN].name = strdup("MPI_Request_free");
  mpi_names[MPI_RECV_INIT_IN].name = strdup("MPI_Recv_init");
  mpi_names[MPI_RECV_INIT_OUT].name = strdup("MPI_Recv_init_out");
  mpi_names[MPI_SEND_INIT_IN].name = strdup("MPI_Send_init");
  mpi_names[MPI_SEND_INIT_OUT].name = strdup("MPI_Send_init_out");
  mpi_names[MPI_GET_ELEMENTS_IN].name = strdup("MPI_Get_elements");
  mpi_names[MPI_GET_COUNT_IN].name = strdup("MPI_Get_count");
  mpi_names[MPI_IBSEND_IN].name = strdup("MPI_Ibsend");
  mpi_names[MPI_IPROBE_IN].name = strdup("MPI_Iprobe");
  mpi_names[MPI_IRSEND_IN].name = strdup("MPI_Irsend");
  mpi_names[MPI_ISSEND_IN].name = strdup("MPI_Issend");
  mpi_names[MPI_PACK_IN].name = strdup("MPI_Pack");
  mpi_names[MPI_PACK_SIZE_IN].name = strdup("MPI_Pack_size");
  mpi_names[MPI_PROBE_IN].name = strdup("MPI_Probe");
  mpi_names[MPI_RSEND_IN].name = strdup("MPI_Rsend");
  mpi_names[MPI_RSEND_INIT_IN].name = strdup("MPI_Rsend_init");
  mpi_names[MPI_RSEND_INIT_OUT].name = strdup("MPI_Rsend_init_out");
  mpi_names[MPI_SENDRECV_IN].name = strdup("MPI_Sendrecv");
  mpi_names[MPI_SENDRECV_REPLACE_IN].name = strdup("MPI_Sendrecv_replace");
  mpi_names[MPI_SSEND_IN].name = strdup("MPI_Ssend");
  mpi_names[MPI_SSEND_INIT_IN].name = strdup("MPI_Ssend_init");
  mpi_names[MPI_SSEND_INIT_OUT].name = strdup("MPI_Ssend_init_out");
  mpi_names[MPI_START_IN].name = strdup("MPI_Start");
  mpi_names[MPI_STARTALL_IN].name = strdup("MPI_Startall");
  mpi_names[MPI_TEST_IN].name = strdup("MPI_Test");
  mpi_names[MPI_TESTALL_IN].name = strdup("MPI_Testall");
  mpi_names[MPI_TESTANY_IN].name = strdup("MPI_Testany");
  mpi_names[MPI_TEST_CANCELLED_IN].name = strdup("MPI_Test_cancelled");
  mpi_names[MPI_TESTSOME_IN].name = strdup("MPI_Testsome");
  mpi_names[MPI_TYPE_COMMIT_IN].name = strdup("MPI_Type_commit");
  mpi_names[MPI_TYPE_CONTIGUOUS_IN].name = strdup("MPI_Type_contiguous");
  mpi_names[MPI_TYPE_EXTENT_IN].name = strdup("MPI_Type_extent");
  mpi_names[MPI_TYPE_FREE_IN].name = strdup("MPI_Type_free");
  mpi_names[MPI_TYPE_HINDEXED_IN].name = strdup("MPI_Type_hindexed");
  mpi_names[MPI_TYPE_HVECTOR_IN].name = strdup("MPI_Type_hvector");
  mpi_names[MPI_TYPE_INDEXED_IN].name = strdup("MPI_Type_indexed");
  mpi_names[MPI_TYPE_INDEXED_OUT].name = strdup("MPI_Type_indexed_out");
  mpi_names[MPI_TYPE_LB_IN].name = strdup("MPI_Type_lb");
  mpi_names[MPI_TYPE_SIZE_IN].name = strdup("MPI_Type_size");
  mpi_names[MPI_TYPE_STRUCT_IN].name = strdup("MPI_Type_struct");
  mpi_names[MPI_TYPE_UB_IN].name = strdup("MPI_Type_ub");
  mpi_names[MPI_TYPE_VECTOR_IN].name = strdup("MPI_Type_vector");
  mpi_names[MPI_UNPACK_IN].name = strdup("MPI_Unpack");
  mpi_names[MPI_WAITALL_IN].name = strdup("MPI_Waitall");
  mpi_names[MPI_WAITANY_IN].name = strdup("MPI_Waitany");
  mpi_names[MPI_WAITSOME_IN].name = strdup("MPI_Waitsome");
  mpi_names[MPI_CART_COORDS_IN].name = strdup("MPI_Cart_coords");
  mpi_names[MPI_CART_CREATE_IN].name = strdup("MPI_Cart_create");
  mpi_names[MPI_CART_GET_IN].name = strdup("MPI_Cart_get");
  mpi_names[MPI_CART_MAP_IN].name = strdup("MPI_Cart_map");
  mpi_names[MPI_CART_SHIFT_IN].name = strdup("MPI_Cart_shift");
  mpi_names[MPI_CARTDIM_GET_IN].name = strdup("MPI_Cartdim_get");
  mpi_names[MPI_DIMS_CREATE_IN].name = strdup("MPI_Dims_create");
  mpi_names[MPI_GRAPH_CREATE_IN].name = strdup("MPI_Graph_create");
  mpi_names[MPI_GRAPH_GET_IN].name = strdup("MPI_Graph_get");
  mpi_names[MPI_GRAPH_MAP_IN].name = strdup("MPI_Graph_map");
  mpi_names[MPI_GRAPH_NEIGHBORS_IN].name = strdup("MPI_Graph_neighbors");
  mpi_names[MPI_GRAPH_NEIGHBORS_COUNT_IN].name = strdup("MPI_Graph_neighbors_count");
  mpi_names[MPI_GRAPHDIMS_GET_IN].name = strdup("MPI_Graphdims_get");
  mpi_names[MPI_TOPO_TEST_IN].name = strdup("MPI_Topo_test");
  mpi_names[MPI_RECV_IDLE_IN].name = strdup("MPI_Recv_idle");
  mpi_names[MPI_CART_RANK_IN].name = strdup("MPI_Cart_rank");
  mpi_names[MPI_CART_SUB_IN].name = strdup("MPI_Cart_sub");
  mpi_names[MPI_FINALIZE_IN].name = strdup("MPI_Finalize");
#endif
  
  for (i = 0; mpi_evt_names[i].name != NULL; i++) {
    mpi_names[mpi_evt_names[i].id].name = strdup(mpi_evt_names[i].name);
    mpi_names[mpi_evt_names[i].id].data = &(mpi_evt_names[i]);
  }

  for (i = 0; olam_evt_names[i].name != NULL; i++) {
    mpi_names[olam_evt_names[i].id].name = strdup(olam_evt_names[i].name);
    mpi_names[olam_evt_names[i].id].data = &(olam_evt_names[i]);
  }
  
  for (i = 0; pvfs_evt_names[i].name != NULL; i++) {
    mpi_names[pvfs_evt_names[i].id].name = strdup(pvfs_evt_names[i].name);
    mpi_names[pvfs_evt_names[i].id].data = &(pvfs_evt_names[i]);
  }
}

char *name_get(int id)
{
  return mpi_names[id].name;
}

evt_name_t *data_get(int id) {
  return (evt_name_t *)(mpi_names[id].data);
}


#ifdef UNIT_TEST
#warning "Unit Testing!!!"

#ifdef OLAM_TEST
#warning "WORKS"
#endif

int main(int argc, char ** argv)
{
  int i;

  name_init();

  for (i = 0; i < MAX_AKY_STATE_NAMES; i++) {
    char *name = name_get(i);
    if (name) {
      evt_name_t *data;
      data = data_get(i);
      
      fprintf(stdout,"Name %s has id %d\n",name,i);
      
      //if ((OLAM_EVT_BASE <= i) && (i <= OLAM_EVT_MAX)) {
      if (data){
        fprintf(stdout,"\tshort: %s\n\tstart_name: %s\n\tstart_id: %d\n\tevt_type:%d\n",
            data->short_name,
            data->start_name,
            data->start_id,
            data->type);
      }
    }
  }
}

#endif
