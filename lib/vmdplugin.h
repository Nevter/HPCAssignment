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
  *      $RCSfile: vmdplugin.h,v $
  *      $Author: johns $       $Locker:  $             $State: Exp $
  *      $Revision: 1.34 $       $Date: 2018/05/02 03:12:56 $
  *
  ***************************************************************************/
 
 #ifndef VMD_PLUGIN_H
 #define VMD_PLUGIN_H
 
 
 /* 
  * Preprocessor tricks to make it easier for us to redefine the names of
  * functions when building static plugins.
  */
 #if !defined(VMDPLUGIN)
 
 #define VMDPLUGIN vmdplugin
 #endif
 
 #define xcat(x, y) cat(x, y)
 
 #define cat(x, y) x ## y 
 
 /*
  *  macros to correctly define plugin function names depending on whether 
  *  the plugin is being compiled for static linkage or dynamic loading. 
  *  When compiled for static linkage, each plugin needs to have unique
  *  function names for all of its entry points.  When compiled for dynamic
  *  loading, the plugins must name their entry points consistently so that
  *  the plugin loading mechanism can find the register, register_tcl, init,
  *  and fini routines via dlopen() or similar operating system interfaces.
  */
 #define VMDPLUGIN_register     xcat(VMDPLUGIN, _register)
 #define VMDPLUGIN_register_tcl xcat(VMDPLUGIN, _register_tcl)
 #define VMDPLUGIN_init         xcat(VMDPLUGIN, _init)
 #define VMDPLUGIN_fini         xcat(VMDPLUGIN, _fini)
 
 #if (defined(WIN32)) 
 #define WIN32_LEAN_AND_MEAN
 #include <windows.h>
 
 #define STATIC_PLUGIN
 
 #if !defined(STATIC_PLUGIN)
 #if defined(VMDPLUGIN_EXPORTS)
 
 BOOL APIENTRY DllMain( HANDLE hModule,
                        DWORD ul_reason_for_call,
                        LPVOID lpReserved
                      )
 {
   return TRUE;
 }
 
 #define VMDPLUGIN_API __declspec(dllexport)
 #else
 #define VMDPLUGIN_API __declspec(dllimport)
 #endif /* VMDPLUGIN_EXPORTS */
 #else  /* ! STATIC_PLUGIN */
 #define VMDPLUGIN_API
 #endif /* ! STATIC_PLUGIN */
 #else
 
 #define VMDPLUGIN_API 
 #endif
 
 #ifdef __cplusplus
 #define VMDPLUGIN_EXTERN extern "C" VMDPLUGIN_API
 #else
 #define VMDPLUGIN_EXTERN extern VMDPLUGIN_API
 #endif  /* __cplusplus */
 
 /* 
  * Plugin API functions start here 
  */
 
 
 VMDPLUGIN_EXTERN int VMDPLUGIN_init(void);
 
 #define vmdplugin_HEAD \
   int abiversion; \
   const char *type; \
   const char *name; \
   const char *prettyname; \
   const char *author; \
   int majorv; \
   int minorv; \
   int is_reentrant; 
 
 typedef struct {
   vmdplugin_HEAD
 } vmdplugin_t;
 
 #define vmdplugin_ABIVERSION  18
 
 #define VMDPLUGIN_THREADUNSAFE 0
 #define VMDPLUGIN_THREADSAFE   1
 
 #define VMDPLUGIN_SUCCESS      0
 #define VMDPLUGIN_ERROR       -1
 
 typedef int (*vmdplugin_register_cb)(void *, vmdplugin_t *);
 
 VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *, vmdplugin_register_cb);
 
 VMDPLUGIN_EXTERN int VMDPLUGIN_register_tcl(void *, void *tcl_interp, 
     vmdplugin_register_cb);
 
 VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void);
 
 #endif   /* VMD_PLUGIN_H */