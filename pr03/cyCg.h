// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
///
/// \file		cyCg.h 
/// \author		Cem Yuksel
/// \version	1.3
/// \date		September 5, 2007
///
/// \brief Cg classes for OpenGL
///
///
/// Add the following line into a cpp file.
///
/// cyCgContext CG;
///
/// Initialize Cg after creating OpenGL context by calling cyCgCreate function.
///
//-------------------------------------------------------------------------------

#ifndef _CY_CG_H_INCLUDED_
#define _CY_CG_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cg/cg.h>
#include <cg/cggl.h>

//-------------------------------------------------------------------------------

/// Keeps the Cg context and profiles.
struct cyCgContext
{
	CGcontext	context;
	CGprofile	vertexProfile, fragmentProfile;
};

/// Add the following line into a source file.
///
/// cyCgContext CG;
extern cyCgContext CG;

//-------------------------------------------------------------------------------

/// When there is an error, this function will print out an error message and exit
inline void cyCgErrorCallback()
{
	CGerror err = cgGetError();
	if (err != CG_NO_ERROR)
	{
		fprintf(stdout,"\n%s\n", cgGetErrorString(err));
		if ( cgIsContext( CG.context ) ) {
			const char *str = cgGetLastListing( CG.context );
			if ( str ) fprintf(stdout,"%s\n",str);
			cgDestroyContext(CG.context);
		}
		exit(0);
	}
}

/// Call this function once after you create the OpenGL context to initialize Cg.
inline void cyCgCreate()
{
	if ( ! cgIsContext( CG.context ) ) {
		cgSetErrorCallback(cyCgErrorCallback);

		// Create context
		CG.context = cgCreateContext();

		// Initialize profiles and compiler options
		CG.vertexProfile = cgGLGetLatestProfile( CG_GL_VERTEX );
		cgGLSetOptimalOptions( CG.vertexProfile );
		CG.fragmentProfile = cgGLGetLatestProfile( CG_GL_FRAGMENT );
		cgGLSetOptimalOptions( CG.fragmentProfile );
	}
}

inline void cyCgEnableVP() { cgGLEnableProfile( CG.vertexProfile ); }	///< Enables Cg vertex program
inline void cyCgEnableFP() { cgGLEnableProfile( CG.fragmentProfile ); } ///< Enables Cg fragment program
inline void cyCgEnable() { cyCgEnableVP(); cyCgEnableFP(); }			///< Enables both vertex and fragment programs

inline void cyCgDisableVP() { cgGLDisableProfile( CG.vertexProfile ); }		///< Disables Cg vertex program
inline void cyCgDisableFP() { cgGLDisableProfile( CG.fragmentProfile ); }	///< Disables Cg fragment program
inline void cyCgDisable() { cyCgDisableVP(); cyCgDisableFP(); }				///< Disables both vertex and fragment programs

inline CGcontext cyCgGetContext() { return CG.context; }					///< Returns the Cg context
inline CGprofile cyCgGetVertexProfile() { return CG.vertexProfile; }		///< Returns the Cg vertex profile
inline CGprofile cyCgGetFragmentProfile() { return CG.fragmentProfile; }	///< Returns the Cg fragment profile


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/// Base class for Cg programs classes.
/// For basic shaders with few uniform parameters, you can use this class directly.
/// If your shader has a lot of uniform parameters, you may want to derive a class
/// from this class, and use OnCreate and OnBind events to initialize the custom
/// parameters. See cyCodeBase solutions section for an example.

class cyCgBaseProgram
{
public:
	cyCgBaseProgram() { program = NULL; params = NULL; progType=0; }
	virtual ~cyCgBaseProgram() { Destroy();  }
	
	//////////////////////////////////////////////////////////////////////////
	///@name Create, destroy and bind methods

	/// This method creates the program from a file.
	/// The argument type can be CG_GL_VERTEX or CG_GL_FRAGMENT, and it determines
	/// if this program is a vertex program or a fragment program.
	/// The argument program_file is the file name of the Cg shader,
	/// and the argument entry is the name of the shader program.
	/// If you would like to use compiler options, you can send a string using args argument,
	/// such that each option must be separated by '|' symbol. For example, the args string
	/// "-DLIMIT=0.02|-v" defines the macro LIMIT which is equal to 0.02 and -v option
	/// prints the computer's version to stdout.
	void CreateFromFile( int type, const char *program_file, const char *entry, const char *args=NULL )
	{
		printf("Compiling %s in %s...", entry, program_file );
		fflush(stdout);
		if ( type == CG_GL_VERTEX ) { CreateFromFile(program_file,cyCgGetVertexProfile(),entry,args); progType=type; }
		if ( type == CG_GL_FRAGMENT ) { CreateFromFile(program_file,cyCgGetFragmentProfile(),entry,args); progType=type; } 
		printf("done.\n");
	}

	/// Destroys created program and deletes registered parameters.
	void Destroy() { DestroyProgram(); DeleteRegisteredParams(); }

	/// Returns true if the program is created and not destroyed.
	CGbool IsCreated() { return cgIsProgram(program); }

	/// Call this function to bind the Cg program. When Cg is enabled, the binded program will execute.
	void Bind() { cgGLBindProgram( program ); OnBind(); }

	CGprogram	GetCgProgram() { return program; }	///< Returns the CGprogram
	int			GetProgType() { return progType; }	///< returns program type (0 = no program defined, CG_GL_VERTEX, or CG_GL_FRAGMENT)



	//////////////////////////////////////////////////////////////////////////
	///@name Accessing uniform parameters by name

	/// Assign value to the uniform parameter with the given name.
	/// These functions can be slow. If you need to change the parameter frequently,
	/// first register the parameter using RegisterParameters, and then use SetRegisteredParameter methods.
	void SetNamedParameter1f( char *name, float x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter1f( param, x ); }
	void SetNamedParameter2f( char *name, float x, float y )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter2f( param, x, y ); }
	void SetNamedParameter3f( char *name, float x, float y, float z )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter3f( param, x, y, z ); }
	void SetNamedParameter4f( char *name, float x, float y, float z, float w )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter4f( param, x, y, z, w ); }
	void SetNamedParameter1fv( char *name, const float *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter1fv( param, x ); }
	void SetNamedParameter2fv( char *name, const float *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter2fv( param, x ); }
	void SetNamedParameter3fv( char *name, const float *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter3fv( param, x ); }
	void SetNamedParameter4fv( char *name, const float *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter4fv( param, x ); }
	void SetNamedParameter1d( char *name, double x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter1d( param, x ); }
	void SetNamedParameter2d( char *name, double x, double y )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter2d( param, x, y ); }
	void SetNamedParameter3d( char *name, double x, double y, double z )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter3d( param, x, y, z ); }
	void SetNamedParameter4d( char *name, double x, double y, double z, double w )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter4d( param, x, y, z, w ); }
	void SetNamedParameter1dv( char *name, const double *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter1dv( param, x ); }
	void SetNamedParameter2dv( char *name, const double *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter2dv( param, x ); }
	void SetNamedParameter3dv( char *name, const double *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter3dv( param, x ); }
	void SetNamedParameter4dv( char *name, const double *x )
			{ CGparameter param = cgGetNamedParameter( program, name ); cgGLSetParameter4dv( param, x ); }

	/// Returns the parameter with the given name
	CGparameter GetNamedParam( char *name ) { return cgGetNamedParameter( program, name ); }



	//////////////////////////////////////////////////////////////////////////
	///@name Accessing uniform parameters by id

	/// Use this function to register the parameters that will be updated frequently.
	/// If you will not update the parameters frequently, you may use SetNamedParameter methods.
	/// Usage: RegisterParameters( 3, "scale|zoom|matrix" );
	/// or
	/// RegisterParameters( 3, "scale zoom matrix" );
	/// This overwrites previous registered parameters (if any).
	void RegisterParameters( int num, char *names )
	{
		DeleteRegisteredParams();
		numParam = num;
		params = new CGparameter[num];
		char *name = strtok (names,"| ");
		int i=0;
		while (name != NULL && i<num) {
			if ( strlen(name) > 0 ) {
				params[i++] = cgGetNamedParameter( program, name );
			}
			if ( i >= numParam ) break;
			name = strtok (NULL, "| ");
		}
	}

	/// Deletes the registered parameters
	void DeleteRegisteredParams()
	{
		if(params) delete [] params;
		params=NULL;
		numParam = 0;
	}

	/// Returns the number of registered parameters
	int GetNumRegisteredParameters() { return numParam; }

	/// Assigns value to the uniform parameter specified by id.
	/// The parameter should be registered first using RegisterParameters method.
	/// The argument id should be smaller than the number of registered parameters.
	void SetRegisteredParameter1f( unsigned int id, float x ) { cgGLSetParameter1f( params[id], x ); }
	void SetRegisteredParameter2f( unsigned int id, float x, float y ) { cgGLSetParameter2f( params[id], x, y ); }
	void SetRegisteredParameter3f( unsigned int id, float x, float y, float z ) { cgGLSetParameter3f( params[id], x, y, z ); }
	void SetRegisteredParameter4f( unsigned int id, float x, float y, float z, float w ) { cgGLSetParameter4f( params[id], x, y, z, w ); }
	void SetRegisteredParameter1fv( unsigned int id, const float *x ) { cgGLSetParameter1fv( params[id], x ); }
	void SetRegisteredParameter2fv( unsigned int id, const float *x ) { cgGLSetParameter2fv( params[id], x ); }
	void SetRegisteredParameter3fv( unsigned int id, const float *x ) { cgGLSetParameter3fv( params[id], x ); }
	void SetRegisteredParameter4fv( unsigned int id, const float *x ) { cgGLSetParameter4fv( params[id], x ); }
	void SetRegisteredParameter1d( unsigned int id, double x ) { cgGLSetParameter1d( params[id], x ); }
	void SetRegisteredParameter2d( unsigned int id, double x, double y ) { cgGLSetParameter2d( params[id], x, y ); }
	void SetRegisteredParameter3d( unsigned int id, double x, double y, double z ) { cgGLSetParameter3d( params[id], x, y, z ); }
	void SetRegisteredParameter4d( unsigned int id, double x, double y, double z, double w ) { cgGLSetParameter4d( params[id], x, y, z, w ); }
	void SetRegisteredParameter1dv( unsigned int id, const double *x ) { cgGLSetParameter1dv( params[id], x ); }
	void SetRegisteredParameter2dv( unsigned int id, const double *x ) { cgGLSetParameter2dv( params[id], x ); }
	void SetRegisteredParameter3dv( unsigned int id, const double *x ) { cgGLSetParameter3dv( params[id], x ); }
	void SetRegisteredParameter4dv( unsigned int id, const double *x ) { cgGLSetParameter4dv( params[id], x ); }


protected:

	//////////////////////////////////////////////////////////////////////////
	///@name Events

	/// This method is called when the program is created.
	/// Overload this method to register custom parameters after creating the program.
	/// See cyCodeBase solutions for an example.
	virtual void OnCreate() {}

	/// This method is called whenever the program is binded.
	/// Overload this method to setup the program each time before you use it.
	/// See cyCodeBase solutions for an example.
	virtual void OnBind() {}


	//////////////////////////////////////////////////////////////////////////
	///@name Protected variables

	CGprogram	program;	///< the CGprogram
	CGparameter *params;	///< registered parameters
	int numParam;			///< number of registered parameters

	int progType;			///< program type (0 = no program defined, CG_GL_VERTEX, or CG_GL_FRAGMENT)


private:
	/// \internal
	// Destroys the cg program. This method is called by the destructor.
	void DestroyProgram() { if ( cgIsProgram(program) ) cgDestroyProgram(program); }

	// Creates the Cg program. This method is called the public CreateFromFile method.
	void CreateFromFile( const char *program_file, CGprofile profile, const char *entry, const char *args )
	{
		DestroyProgram();
		// Prepare command arguments
		int argCount = 0;
		char **argArray = NULL;
		if ( args ) {
			size_t argslen = strlen(args);
			char *tempArgs = new char[argslen+1];
			strcpy(tempArgs,args);
			char *arg = strtok (tempArgs,"|");
			int argCount=0;
			size_t maxlen=0;
			while (arg != NULL) {
				size_t len = strlen(arg);
				if ( len > maxlen ) maxlen = len;
				argCount++;
				arg = strtok (NULL, "|");
			}
			argArray = new char*[argCount+1];
			for ( int i=0; i<argCount; i++ ) argArray[i] = new char[maxlen+1];
			argArray[argCount] = NULL;
			strcpy(tempArgs,args);
			arg = strtok (tempArgs,"|");
			int n=0;
			while (arg != NULL) {
				strcpy( argArray[n++], arg );
				arg = strtok (NULL, "|");
			}
			strcpy(tempArgs,args);
			delete [] tempArgs;
		}

		// Create the vertex program
		program = cgCreateProgramFromFile( cyCgGetContext(), CG_SOURCE, program_file, profile, entry, (const char**) argArray);

		if ( args ) {
			for ( int i=0; i<argCount; i++) delete [] argArray[i];
			delete [] argArray;
		}

		// Load the program
		cgGLLoadProgram(program);
		OnCreate();
	}

};

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

/// Multiple versions of the same Cg program.
/// Sometimes you may want to compile several versions of the same Cg program
/// with different compiler options. This template class makes this process easier.
/// It keeps an array of Cg programs. You can compile them all at once,
/// and access them individually.

template <class T> class cyCgProgramContainer
{
public:
	cyCgProgramContainer() { progs=NULL; }
	virtual ~cyCgProgramContainer() { if (progs) delete [] progs; }

	/// Creates a number of programs using the same file and same program entry
	/// with different compiler options. The argument numProgs defines the number of programs,
	/// type is either CG_GL_VERTEX or CG_GL_FRAGMENT, program_file is the file name of the
	/// Cg program, and entry is the name of the shader. This method expects a separate
	/// argument string for each program.
	/// Ex:
	/// CreateFromFile( 2, CG_GL_VERTEX, file, entry, "-DLEVEL=1", "-DLEVEL=2" );
	void CreateFromFile( int numProgs, int type, const char *program_file, const char *entry, ... )
	{
		if ( progs ) delete [] progs;
		progs = new T [numProgs];
		va_list args;
		va_start(args,entry);
		for ( int i=0; i<numProgs; i++ ) {
			const char *arg = va_arg( args, const char* );
			progs[i].CreateFromFile( type, program_file, entry, arg );
		}
		va_end( args );
	}

	/// You can access each individual program using this operator.
	T&		operator [] ( int id ) { return progs[id]; }

	/// Binds the specified Cg program.
	void	Bind( int id ) { progs[id].Bind(); }

protected:
	T	*progs;	// The list of Cg programs
};

/// List of cyCgBaseProgram
typedef cyCgProgramContainer<cyCgBaseProgram> cyCgBaseProgramContainer;

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

namespace cy {
	typedef cyCgContext CgContext;
	typedef cyCgBaseProgram CgBaseProgram;
	template <class T> class CgProgramContainer : public cyCgProgramContainer<T> {};
	typedef cyCgBaseProgramContainer CgBaseProgramContainer;

	static void (*CgCreate)() = cyCgCreate;
	static void (*CgEnableVP)() = cyCgEnableVP;
	static void (*CgEnableFP)() = cyCgEnableFP;
	static void (*CgEnable)() = cyCgEnable;
	static void (*CgDisableVP)() = cyCgDisableVP;
	static void (*CgDisableFP)() = cyCgDisableFP;
	static void (*CgDisable)() = cyCgDisable;
	static CGcontext (*CgGetContext)() = cyCgGetContext;
	static CGprofile (*CgGetVertexProfile)() = cyCgGetVertexProfile;
	static CGprofile (*CgGetFragmentProfile)() = cyCgGetFragmentProfile;
}

//-------------------------------------------------------------------------------

#endif

