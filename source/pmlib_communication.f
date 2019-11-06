!---------------------------------------------------------------------------------!
! pmlib_communication.f
!   version: 3D multi-core
!---------------------------------------------------------------------------------!
MODULE pmlib_mod_communication

USE pmlib_mod_parameters
USE pmlib_mod_write 

IMPLICIT NONE

  TYPE class_comm_buffer
    INTEGER                                 :: nsend, nrecieve
    INTEGER,DIMENSION(pmlib_ndim)           :: ncell_send, ncell_recieve
    INTEGER                                 :: npart_send, npart_recieve

    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: mesh_send    => NULL()
    REAL(MK),DIMENSION(:,:,:,:),POINTER     :: mesh_recieve => NULL()

    REAL(MK),DIMENSION(:,:),POINTER         :: part_send    => NULL()
    REAL(MK),DIMENSION(:,:),POINTER         :: part_recieve => NULL()
  END TYPE class_comm_buffer

  TYPE class_communication
    INTEGER                                 :: ncomm
    INTEGER                                 :: npacked
    INTEGER                                 :: ctype
    INTEGER,DIMENSION(:),POINTER            :: proc_send,proc_recieve
    INTEGER,DIMENSION(:,:),POINTER          :: imin_send, imax_send
    INTEGER,DIMENSION(:,:),POINTER          :: imin_recieve, imax_recieve
    INTEGER,DIMENSION(:,:),POINTER          :: ncell_send, ncell_recieve
    INTEGER,DIMENSION(:),POINTER            :: npart_send, npart_recieve
    INTEGER,DIMENSION(:),POINTER            :: ipart0, ipart1, ipart2
    INTEGER,DIMENSION(:),POINTER            :: comm_dir
    INTEGER                                 :: ncommseq
    INTEGER,DIMENSION(:,:,:),POINTER        :: icommseq
  END TYPE class_communication

  TYPE(class_communication) :: communication

  TYPE(class_comm_buffer),DIMENSION(:),POINTER :: comm_buffer => NULL()

  INTERFACE pmlib_comm_pack
    MODULE PROCEDURE pmlib_comm_pack_mesh_logical
    MODULE PROCEDURE pmlib_comm_pack_mesh_real
    MODULE PROCEDURE pmlib_comm_pack_mesh_complex
    MODULE PROCEDURE pmlib_comm_pack_part
  END INTERFACE

  INTERFACE pmlib_comm_unpack
    MODULE PROCEDURE pmlib_comm_unpack_mesh_logical
    MODULE PROCEDURE pmlib_comm_unpack_mesh_real
    MODULE PROCEDURE pmlib_comm_unpack_mesh_complex
    MODULE PROCEDURE pmlib_comm_unpack_part
    MODULE PROCEDURE pmlib_comm_unpack_part_pos_check
  END INTERFACE

CONTAINS

#include "communication/pmlib_comm_pack.f"

#include "communication/pmlib_comm_send.f"

#include "communication/pmlib_comm_unpack.f"

#include "communication/pmlib_comm_finalise.f"


END MODULE pmlib_mod_communication
