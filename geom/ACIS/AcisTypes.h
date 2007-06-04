

#if CUBIT_ACIS_VERSION < 1000
#define SPAvector       vector
#define SPAunit_vector  unit_vector
#define SPAposition     position
#define SPAmatrix       matrix
#define SPAtransf       transf
#define SPAinterval     interval
#define SPApar_pos      par_pos
#define SPApar_vec      par_vec
#define SPApar_dir      par_dir
#define SPApar_box      par_box
#define SPApar_transf   par_transf
#define SPAbox          box
#define SPAnvector      nvector
#define SPAparameter    parameter
#define SPAresabs       resabs
#define SPAresnor       resnor
#define SPAresfit       resfit
#define SPAresmch       resmch
#endif

#define IS_ENTITY_TYPE(E,T) (E->identity(T##_LEVEL) == T##_TYPE)
