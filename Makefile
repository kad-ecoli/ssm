CXX=g++
CXXFLAGS=-O3 -I${PWD} -I${PWD}/mmdb
LDFLAGS=-static
SSM_HEADERS=ss_csia.h ss_graph.h ssm_align.h ssm_superpose.h ss_vxedge.h
SSM_CPP=ss_csia.cpp ss_graph.cpp ssm_align.cpp ssm_superpose.cpp ss_vxedge.cpp
MMDB_CPP=mmdb/mmdb_manager.cpp  \
         mmdb/mmdb_coormngr.cpp \
	 mmdb/mmdb_selmngr.cpp  \
	 mmdb/mmdb_graph.cpp    \
	 mmdb/mmdb_file.cpp     \
	 mmdb/mmdb_cryst.cpp    \
	 mmdb/mmdb_uddata.cpp   \
	 mmdb/mmdb_mask.cpp     \
	 mmdb/mmdb_symop.cpp    \
	 mmdb/mmdb_utils.cpp    \
	 mmdb/mmdb_title.cpp    \
	 mmdb/mmdb_mmcif.cpp    \
	 mmdb/mmdb_cifdefs.cpp  \
	 mmdb/mmdb_model.cpp    \
	 mmdb/mmdb_chain.cpp    \
	 mmdb/mmdb_bondmngr.cpp \
	 mmdb/mmdb_atom.cpp     \
	 mmdb/mmdb_tables.cpp   \
	 mmdb/bfgs_min.cpp      \
	 mmdb/stream_.cpp       \
	 mmdb/mattype_.cpp      \
	 mmdb/hybrid_36.cpp     \
	 mmdb/linalg_.cpp       \
	 mmdb/math_.cpp         \
	 mmdb/file_.cpp
PROG=ssm

all: ${PROG}

libmmdb.a: ${MMDB_CPP}
	cd mmdb;${CXX} ${CXXFLAGS} -c *.cpp;ar rvs ../libmmdb.a *.o;cd ..

ssm: libmmdb.a ${SSM_HEADERS} ${SSM_CPP} ${MMDB_CPP} ssm.cpp
	${CXX} ${CXXFLAGS} ${SSM_CPP} $@.cpp libmmdb.a -o $@ ${LDFLAGS}

clean-obj:
	rm mmdb/*.o

clean: clean-obj
	rm ${PROG}
