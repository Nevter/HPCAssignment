 /***************************************************************************
  *cr
  *cr            (C) Copyright 1995-2006 The Board of Trustees of the
  *cr                        University of Illinois
  *cr                         All Rights Reserved
  *cr
  ***************************************************************************/
 
 /***************************************************************************
  * RCS INFORMATION:
  *
  *      $RCSfile: molfile_plugin.h,v $
  *      $Author: johns $       $Locker:  $             $State: Exp $
  *      $Revision: 1.110 $       $Date: 2018/05/02 03:12:56 $
  *
  ***************************************************************************/
 
 #ifndef MOL_FILE_PLUGIN_H
 #define MOL_FILE_PLUGIN_H
 
 #include "vmdplugin.h"
 
 #if defined(DESRES_READ_TIMESTEP2)
 /* includes needed for large integer types used for frame counts */
 #include <sys/types.h>
 typedef ssize_t molfile_ssize_t;      
 #endif
 
 #define MOLFILE_PLUGIN_TYPE "mol file reader"
 
 #define MOLFILE_CONVERTER_PLUGIN_TYPE "mol file converter"
 
 /* File plugin symbolic constants for better code readability */
 #define MOLFILE_SUCCESS           0   
 #define MOLFILE_EOF              -1   
 #define MOLFILE_ERROR            -1   
 #define MOLFILE_NOSTRUCTUREDATA  -2   
 #define MOLFILE_NUMATOMS_UNKNOWN -1   
 #define MOLFILE_NUMATOMS_NONE     0   
 #define MOLFILE_BUFSIZ           81   
 #define MOLFILE_BIGBUFSIZ        4096 
 #define MOLFILE_MAXWAVEPERTS     25   
 #define MOLFILE_DIRECTIO_MIN_BLOCK_SIZE 4096
 #define MOLFILE_DIRECTIO_MAX_BLOCK_SIZE 4096
 
 
 typedef struct {
   char database[81];   
   char accession[81];  
   char date[81];       
   char title[81];      
   int remarklen;       
   char *remarks;       
 } molfile_metadata_t;
 
 
 /* 
  * Struct for specifying atoms in a molecular structure.  The first 
  * six components are required, the rest are optional and their presence is 
  * indicating by setting the corresponding bit in optsflag.  When omitted,
  * the application (for read_structure) or plugin (for write_structure) 
  * must be able to supply default values if the missing parameters are 
  * part of its internal data structure.
  * Note that it is not possible to specify coordinates with this structure.
  * This is intentional; all coordinate I/O is done with the read_timestep and 
  * write_timestep functions. 
  */
 
 typedef struct {
   /* these fields absolutely must be set or initialized to empty */
   char name[16];      
   char type[16];      
   char resname[8];    
   int resid;          
   char segid[8];      
 #if 0 && vmdplugin_ABIVERSION > 10000
   /* The new PDB file formats allows for much larger structures, */
   /* which can therefore require longer chain ID strings.  The   */
   /* new PDBx/mmCIF file formats do not have length limits on    */
   /* fields, so PDB chains could be arbitrarily long strings     */
   /* in such files.  At present, we know we need at least 3-char */
   /* chains for existing PDBx/mmCIF files.                       */
   char chain[4];      
 #else
   char chain[2];      
 #endif
   /* rest are optional; use optflags to specify what's present   */
   char altloc[2];     
   char insertion[2];  
   float occupancy;    
   float bfactor;      
   float mass;         
   float charge;       
   float radius;       
   int atomicnumber;   
 #if 0
   char complex[16];
   char assembly[16];
   int qmregion;
   int qmregionlink;
   int qmlayer;
   int qmlayerlink;
   int qmfrag;
   int qmfraglink;
   string qmecp;
   int qmadapt;
   int qmect;          
   int qmparam;
   int autoparam;
 #endif
 
 #if defined(DESRES_CTNUMBER)
   int ctnumber;       
 #endif
 } molfile_atom_t;
 
 #define MOLFILE_NOOPTIONS     0x0000 
 #define MOLFILE_INSERTION     0x0001 
 #define MOLFILE_OCCUPANCY     0x0002 
 #define MOLFILE_BFACTOR       0x0004 
 #define MOLFILE_MASS          0x0008 
 #define MOLFILE_CHARGE        0x0010 
 #define MOLFILE_RADIUS        0x0020 
 #define MOLFILE_ALTLOC        0x0040 
 #define MOLFILE_ATOMICNUMBER  0x0080 
 #define MOLFILE_BONDSSPECIAL  0x0100 
 #if defined(DESRES_CTNUMBER)
 #define MOLFILE_CTNUMBER      0x0200 
 #endif
 #define MOLFILE_BADOPTIONS    0xFFFFFFFF 
 
 
 #define MOLFILE_QMTS_NOOPTIONS     0x0000 
 #define MOLFILE_QMTS_GRADIENT      0x0001 
 #define MOLFILE_QMTS_SCFITER       0x0002
 
 typedef struct molfile_timestep_metadata {
   unsigned int count;                  
   unsigned int avg_bytes_per_timestep; 
   int has_velocities;                  
 } molfile_timestep_metadata_t;
 
 /*
  * Per-timestep atom coordinates and periodic cell information
  */ 
 typedef struct {
   float *coords;        
   float *velocities;    
   float A, B, C, alpha, beta, gamma; 
   double physical_time; 
 #if defined(DESRES_READ_TIMESTEP2)
   /* HACK to support generic trajectory information */
   double total_energy;
   double potential_energy;
   double kinetic_energy;
   double extended_energy;
   double force_energy;
   double total_pressure;
 #endif
 
 } molfile_timestep_t;
 
 
 typedef struct {
   char dataname[256];   
   float origin[3];      
   /*
    * x/y/z axis:
    * These the three cell sides, providing both direction and length
    * (not unit vectors) for the x, y, and z axes.  In the simplest
    * case, these would be <size,0,0> <0,size,0> and <0,0,size) for 
    * an orthogonal cubic volume set.  For other cell shapes these
    * axes can be oriented non-orthogonally, and the parallelpiped
    * may have different side lengths, not just a cube/rhombus.
    */
   float xaxis[3];       
   float yaxis[3];       
   float zaxis[3];       
   /*
    * x/y/z size: 
    * Number of grid cells along each axis.  This is _not_ the
    * physical size of the box, this is the number of voxels in each
    * direction, independent of the shape of the volume set. 
    */
   int xsize;            
   int ysize;            
   int zsize;            
 #if vmdplugin_ABIVERSION > 16
   int has_scalar;       
   int has_gradient;     
   int has_variance;     
 #endif
   int has_color;        
 } molfile_volumetric_t;
 
 
 #if vmdplugin_ABIVERSION > 16
 
 typedef struct {
   int setidx;           
   float *scalar;        
   float *gradient;      
   float *variance;      
   float *rgb3f;         
   unsigned char *rgb3u; 
 } molfile_volumetric_readwrite_t;
 #endif
 
 
 /**************************************************************
  **************************************************************
  ****                                                      ****
  ****          Data structures for QM files                ****
  ****                                                      ****
  **************************************************************
  **************************************************************/
 
 /* macros for the convergence status of a QM calculation. */
 #define MOLFILE_QMSTATUS_UNKNOWN       -1 /* don't know yet */
 #define MOLFILE_QMSTATUS_OPT_CONV       0 /* optimization converged */
 #define MOLFILE_QMSTATUS_SCF_NOT_CONV   1 /* SCF convergence failed */
 #define MOLFILE_QMSTATUS_OPT_NOT_CONV   2 /* optimization not converged */
 #define MOLFILE_QMSTATUS_FILE_TRUNCATED 3 /* file was truncated */
 
 /* macros describing the SCF method (SCFTYP in GAMESS) */
 #define MOLFILE_SCFTYPE_UNKNOWN -1 /* no info about the method  */
 #define MOLFILE_SCFTYPE_NONE     0 /* calculation didn't make use of SCF */
 #define MOLFILE_SCFTYPE_RHF      1 /* restricted Hartree-Fock   */
 #define MOLFILE_SCFTYPE_UHF      2 /* unrestricted Hartree-Fock */
 #define MOLFILE_SCFTYPE_ROHF     3 /* restricted open-shell Hartree-Fock */
 #define MOLFILE_SCFTYPE_GVB      4 /* generalized valence bond orbitals  */
 #define MOLFILE_SCFTYPE_MCSCF    5 /* multi-configuration SCF   */
 #define MOLFILE_SCFTYPE_FF       6 /* classical force-field based sim.   */
 
 /* macros describing the type of calculation (RUNTYP in GAMESS) */
 #define MOLFILE_RUNTYPE_UNKNOWN    0  /* single point run */
 #define MOLFILE_RUNTYPE_ENERGY     1  /* single point run */
 #define MOLFILE_RUNTYPE_OPTIMIZE   2  /* geometry optimization */
 #define MOLFILE_RUNTYPE_SADPOINT   3  /* saddle point search */
 #define MOLFILE_RUNTYPE_HESSIAN    4  /* Hessian/frequency calculation */
 #define MOLFILE_RUNTYPE_SURFACE    5  /* potential surface scan */
 #define MOLFILE_RUNTYPE_GRADIENT   6  /* energy gradient calculation */
 #define MOLFILE_RUNTYPE_MEX        7  /* minimum energy crossing */
 #define MOLFILE_RUNTYPE_DYNAMICS   8  /* Any type of molecular dynamics
                                        * e.g. Born-Oppenheimer, Car-Parinello,
                                        * or classical MD */
 #define MOLFILE_RUNTYPE_PROPERTIES 9  /* Properties were calculated from a
                                        * wavefunction that was read from file */
 
 
 typedef struct {
   /* hessian data */
   int nimag;                    
   int nintcoords;               
   int ncart;                    
   /* orbital/basisset data */
   int num_basis_funcs;          
   int num_basis_atoms;          
   int num_shells;               
   int wavef_size;               
   /* everything else */
   int have_sysinfo;
   int have_carthessian;         
   int have_inthessian;          
   int have_normalmodes;         
 } molfile_qm_metadata_t;
 
 
 typedef struct {
   int nproc;             
   int memory;            
   int runtype;           
   int scftype;           
   int status;            
   int num_electrons;     
   int totalcharge;       
   int num_occupied_A;    
   int num_occupied_B;    
   double *nuc_charge;    
   char basis_string[MOLFILE_BUFSIZ];    
   char runtitle[MOLFILE_BIGBUFSIZ];     
   char geometry[MOLFILE_BUFSIZ];        
   char version_string[MOLFILE_BUFSIZ];  
 } molfile_qm_sysinfo_t;
 
 
 typedef struct {
   int *num_shells_per_atom; 
   int *num_prim_per_shell;  
   float *basis;             
   int *atomic_number;       
   int *angular_momentum;    
   int *shell_types;         
 } molfile_qm_basis_t;
 
 
 typedef struct {
   double *carthessian;  
   int    *imag_modes;   
   double *inthessian;   
   float *wavenumbers;   
   float *intensities;   
   float *normalmodes;   
 } molfile_qm_hessian_t;
 
 
 typedef struct {
   molfile_qm_sysinfo_t run;             /* system info  */
   molfile_qm_basis_t   basis;           /* basis set info */
   molfile_qm_hessian_t hess;            /* hessian info */
 } molfile_qm_t;
 
 
 
 enum molfile_qm_wavefunc_type {
   MOLFILE_WAVE_CANON,    MOLFILE_WAVE_GEMINAL,
   MOLFILE_WAVE_MCSCFNAT, MOLFILE_WAVE_MCSCFOPT,
   MOLFILE_WAVE_CINATUR,
   MOLFILE_WAVE_PIPEK,  MOLFILE_WAVE_BOYS, MOLFILE_WAVE_RUEDEN,
   MOLFILE_WAVE_NAO,    MOLFILE_WAVE_PNAO, MOLFILE_WAVE_NHO, 
   MOLFILE_WAVE_PNHO,   MOLFILE_WAVE_NBO,  MOLFILE_WAVE_PNBO, 
   MOLFILE_WAVE_PNLMO,  MOLFILE_WAVE_NLMO, MOLFILE_WAVE_MOAO, 
   MOLFILE_WAVE_NATO,   MOLFILE_WAVE_UNKNOWN
 };
 
 
 enum molfile_qm_charge_type {
   MOLFILE_QMCHARGE_UNKNOWN,
   MOLFILE_QMCHARGE_MULLIKEN, MOLFILE_QMCHARGE_LOWDIN,
   MOLFILE_QMCHARGE_ESP, MOLFILE_QMCHARGE_NPA
 };
 
 
 
 typedef struct molfile_qm_timestep_metadata {
   unsigned int count;                  
   unsigned int avg_bytes_per_timestep; 
   int has_gradient;                    
   int num_scfiter;                     
   int num_orbitals_per_wavef[MOLFILE_MAXWAVEPERTS]; 
   int has_orben_per_wavef[MOLFILE_MAXWAVEPERTS]; 
   int has_occup_per_wavef[MOLFILE_MAXWAVEPERTS]; 
   int num_wavef ;                      
   int wavef_size;                      
   int num_charge_sets;                 
 } molfile_qm_timestep_metadata_t;
 
 
 typedef struct {
   int   type;               
   int   spin;               
   int   excitation;         
   int   multiplicity;       
   char info[MOLFILE_BUFSIZ]; 
   double energy;            
   float *wave_coeffs;       
   float *orbital_energies;  
   float *occupancies;       
   int   *orbital_ids;       
 } molfile_qm_wavefunction_t;
 
 
 typedef struct {
   molfile_qm_wavefunction_t *wave; 
   float  *gradient;         
   double *scfenergies;      
   double *charges;          
   int    *charge_types;     
 } molfile_qm_timestep_t;
 
 
 /**************************************************************
  **************************************************************/
 
 
 
 
 enum molfile_graphics_type {
   MOLFILE_POINT,  MOLFILE_TRIANGLE, MOLFILE_TRINORM, MOLFILE_NORMS, 
   MOLFILE_LINE,   MOLFILE_CYLINDER, MOLFILE_CAPCYL,  MOLFILE_CONE,    
   MOLFILE_SPHERE, MOLFILE_TEXT,     MOLFILE_COLOR,   MOLFILE_TRICOLOR
 };
 
 typedef struct {
   int type;             /* One of molfile_graphics_type */
   int style;            /* A general style parameter    */
   float size;           /* A general size parameter     */
   float data[9];        /* All data for the element     */
 } molfile_graphics_t;
 
 
 /*
  * Types for raw graphics elements stored in files.  Data for each type
  * should be stored by the plugin as follows:
 
 type        data                                     style       size
 ----        ----                                     -----       ----
 point       x, y, z                                              pixel size
 triangle    x1,y1,z1,x2,y2,z2,x3,y3,z3                 
 trinorm     x1,y1,z1,x2,y2,z2,x3,y3,z3                 
             the next array element must be NORMS
 tricolor    x1,y1,z1,x2,y2,z2,x3,y3,z3                 
             the next array elements must be NORMS
             the following element must be COLOR, with three RGB triples
 norms       x1,y1,z1,x2,y2,z2,x3,y3,z3                 
 line        x1,y1,z1,x2,y2,z2                        0=solid     pixel width
                                                      1=stippled
 cylinder    x1,y1,z1,x2,y2,z2                        resolution  radius
 capcyl      x1,y1,z1,x2,y2,z2                        resolution  radius
 sphere      x1,y1,z1                                 resolution  radius
 text        x, y, z, up to 24 bytes of text                      pixel size
 color       r, g, b
 */
 
 
 typedef struct {
   vmdplugin_HEAD
 
   const char *filename_extension;
 
   void *(* open_file_read)(const char *filepath, const char *filetype, 
       int *natoms);
   
   int (*read_structure)(void *, int *optflags, molfile_atom_t *atoms);
 
   int (*read_bonds)(void *, int *nbonds, int **from, int **to, float **bondorder, 
                     int **bondtype, int *nbondtypes, char ***bondtypename);
 
   int (* read_next_timestep)(void *, int natoms, molfile_timestep_t *);
 
   void (* close_file_read)(void *);
    
   void *(* open_file_write)(const char *filepath, const char *filetype, 
       int natoms);
   
   int (* write_structure)(void *, int optflags, const molfile_atom_t *atoms);
 
   int (* write_timestep)(void *, const molfile_timestep_t *);
   
   void (* close_file_write)(void *);
 
   int (* read_volumetric_metadata)(void *, int *nsets, 
         molfile_volumetric_t **metadata);
 
   int (* read_volumetric_data)(void *, int set, float *datablock, 
         float *colorblock);
 #if vmdplugin_ABIVERSION > 16
   int (* read_volumetric_data_ex)(void *, molfile_volumetric_readwrite_t *v);
 #endif
 
   int (* read_rawgraphics)(void *, int *nelem, const molfile_graphics_t **data);
 
   int (* read_molecule_metadata)(void *, molfile_metadata_t **metadata);
   
   int (* write_bonds)(void *, int nbonds, int *from, int *to, float *bondorder, 
                      int *bondtype, int nbondtypes, char **bondtypename);
 
   int (* write_volumetric_data)(void *, molfile_volumetric_t *metadata,
                                 float *datablock, float *colorblock);
 #if vmdplugin_ABIVERSION > 16
   int (* write_volumetric_data_ex)(void *, molfile_volumetric_t *metadata,
                                    molfile_volumetric_readwrite_t *v);
 #endif
 
   int (* read_angles)(void *handle, int *numangles, int **angles, int **angletypes,
                       int *numangletypes, char ***angletypenames, int *numdihedrals,
                       int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
                       char ***dihedraltypenames, int *numimpropers, int **impropers,        
                       int **impropertypes, int *numimpropertypes, char ***impropertypenames,
                       int *numcterms, int **cterms, int *ctermcols, int *ctermrows);
 
   int (* write_angles)(void *handle, int numangles, const int *angles, const int *angletypes,
                        int numangletypes, const char **angletypenames, int numdihedrals,
                        const int *dihedrals, const int *dihedraltypes, int numdihedraltypes,
                        const char **dihedraltypenames, int numimpropers, 
                        const int *impropers, const int *impropertypes, int numimpropertypes,
                        const char **impropertypenames, int numcterms,  const int *cterms, 
                        int ctermcols, int ctermrows);
 
 
   int (* read_qm_metadata)(void *, molfile_qm_metadata_t *metadata);
 
 
   int (* read_qm_rundata)(void *, molfile_qm_t *qmdata);
 
 
 #if vmdplugin_ABIVERSION > 17
   int (* read_timestep_pagealign_size)(void *, int *pagealignsize);
 #endif
 
 
   int (* read_timestep)(void *, int natoms, molfile_timestep_t *,
                         molfile_qm_metadata_t *, molfile_qm_timestep_t *);
 
   int (* read_timestep_metadata)(void *, molfile_timestep_metadata_t *);
   int (* read_qm_timestep_metadata)(void *, molfile_qm_timestep_metadata_t *);
 
 #if defined(DESRES_READ_TIMESTEP2)
 
   int (* read_timestep2)(void *, molfile_ssize_t index, molfile_timestep_t *);
 
   molfile_ssize_t (* read_times)( void *,
                                   molfile_ssize_t start,
                                   molfile_ssize_t count,
                                   double * times );
 #endif
 
   int (* cons_fputs)(const int, const char*);
 
 } molfile_plugin_t;
 
 #endif
 