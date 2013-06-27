# cython: profile=True
import random
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport bool

__all__ = [ "Spb", "Chromosome",
           "Centromere", "PlugSite", "Spindle"]

ITYPE = np.int
CTYPE = np.float
ctypedef np.int_t ITYPE_t
ctypedef np.float_t CTYPE_t

cdef class Organite
cdef class Spb
cdef class Chromosome
cdef class Centromere
cdef class PlugSite
cdef class Spindle

cdef class Spindle(object) :
    cdef public object KD
    cdef public object all_plugsites

cdef class Organite(object):
    cdef public int num_steps
    cdef public np.ndarray traj
    cdef public float pos
    cdef public object parent, KD
    cdef void set_pos(Organite, float, int time_point=*)
    cdef float get_pos(Organite, int time_point=*)


cdef class Spb(Organite):
    cdef public int side

cdef class Chromosome(Organite):
    cdef public Centromere cen_A, cen_B
    cdef public np.ndarray correct_history, erroneous_history
    cdef bool is_right_A(Chromosome)
    cdef int delta(Chromosome)
    cdef float pair_dist(Chromosome)
    cdef float plug_dist(Chromosome, float)
    cdef float center(Chromosome)
    cdef np.ndarray center_traj(Chromosome)
    cdef bool at_rightpole(Chromosome, float)
    cdef bool at_leftpole(Chromosome, float)
    cdef bool at_pole(Chromosome, int side=*, float tol=*)
    cdef public int ch_id
    cdef public int id


cdef class Centromere(Organite):
    cdef public str tag
    cdef public Chromosome chromosome
    cdef public float toa
    cdef public np.ndarray plug_vector
    cdef public list plugsites
    cdef void calc_plug_vector(Centromere)
    cdef float P_attachleft(Centromere)
    cdef int left_plugged(Centromere)
    cdef int right_plugged(Centromere)
    cdef bool at_rightpole(Centromere, float tol=*)
    cdef bool at_leftpole(Centromere, float tol=*)

cdef class PlugSite(Organite):
    cdef public Centromere centromere
    cdef public str tag
    cdef public int plug_state, plugged
    cdef public np.ndarray state_hist
    cdef public float P_att
    cdef public void set_plug_state(PlugSite, int, int time_point=*)
    cdef public float calc_ldep(PlugSite)
    cdef public void plug_unplug(PlugSite, int)
    cdef public float P_det(PlugSite)
    cdef public int site_id
