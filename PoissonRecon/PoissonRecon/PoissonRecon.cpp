
#ifdef WIN32
#include <windows.h>
#include <Psapi.h>
#endif
#include <common/interfaces.h>
#include <cstdio>
#include "Src/MyTime.h"
#include "Src/MemoryUsage.h"
#include "Src/MarchingCubes.h"
#include "Src/Octree.h"
#include "Src/SparseMatrix.h"
#include "Src/CmdLineParser.h"
#include "Src/PPolynomial.h"
//#include "Src/Ply.h"
void DumpOutput( const char* format , ... );
void DumpOutput2( char* str , const char* format , ... );

#include "Src/PointStream.h"

#include "Src/MultiGridOctreeData.h"
#include "PoissonRecon.h"
#include "NearestNeighborSearches.h"
#include <QtScript>


#define zoom_mm_distance 0.07
#define min_box_size 30

void DumpOutput( const char* format , ... )
{
  char buf[4096];
  va_list marker;
  va_start( marker, format );

  vsprintf(buf,format,marker);
  va_end( marker );

  qDebug(buf);
 }
void DumpOutput2(std::vector< char* >& comments  , const char* format , ... )
{
  char buf[4096];
  va_list marker;
  va_start( marker, format );

  vsprintf(buf,format,marker);
  va_end( marker );
  qDebug(buf);
}


#if defined( _WIN32 ) || defined( _WIN64 )
double PeakMemoryUsageMB( void )
{
	HANDLE h = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	return GetProcessMemoryInfo( h , &pmc , sizeof(pmc) ) ? ( (double)pmc.PeakWorkingSetSize )/(1<<20) : 0;
}
#endif // _WIN32 || _WIN64

#if defined( _WIN32 ) || defined( _WIN64 )
inline double to_seconds( const FILETIME& ft )
{
	const double low_to_sec=100e-9; // 100 nanoseconds
	const double high_to_sec=low_to_sec*4294967296.0;
	return ft.dwLowDateTime*low_to_sec+ft.dwHighDateTime*high_to_sec;
}
#endif // _WIN32 || _WIN64


template< class Real >
struct OctreeProfiler
{
	Octree< Real >& tree;
	double t;

	OctreeProfiler( Octree< Real >& t ) : tree(t) { ; }
	void start( void ){ t = Time() , tree.resetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}
	void dumpOutput( const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) DumpOutput( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         DumpOutput(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) DumpOutput( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         DumpOutput(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}
	void dumpOutput2( std::vector< char* >& comments , const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) DumpOutput2( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         DumpOutput2( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) DumpOutput2( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         DumpOutput2( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}
};

PoissonReconstruction::PoissonReconstruction()
{
}
template <class Real>
class PoissonParam
{
public:
  int MaxDepthVal;
  int MaxSolveDepthVal;
  int KernelDepthVal;
  int MinDepthVal;
  int FullDepthVal;
  Real SamplesPerNodeVal;
  Real ScaleVal;
  bool ConfidenceFlag;
  bool CleanFlag;
  bool DensityFlag;
  Real PointWeightVal;
  int AdaptiveExponentVal;
  int BoundaryTypeVal;
  bool CompleteFlag;
  bool NonManifoldFlag;
  bool ShowResidualFlag;
  int CGDepthVal;
  int ItersVal;
  Real CSSolverAccuracyVal;
  
  bool VerboseFlag;
  int ThreadsVal;
  bool LinearFitFlag;
  float LowResIterMultiplierVal;
  float ColorVal;
  
  PoissonParam()
  {
    MaxDepthVal=8;
    MaxSolveDepthVal=-1;
    KernelDepthVal=-1;
    MinDepthVal=0;
    FullDepthVal=5;
    SamplesPerNodeVal =1.5f;
    ScaleVal=1.1f;
    ConfidenceFlag=false;
    CleanFlag=false;
    DensityFlag=false;
    PointWeightVal = 4.0f;
    AdaptiveExponentVal=1;
    BoundaryTypeVal=1;
    CompleteFlag=false;
    NonManifoldFlag=false;
    ShowResidualFlag=false;
    CGDepthVal=0;
    ItersVal=8;
    CSSolverAccuracyVal=1e-3f;
    
    VerboseFlag=true;
    ThreadsVal=omp_get_num_procs();
    LinearFitFlag = false;
    LowResIterMultiplierVal=1.f;
    ColorVal=16.0f;
  }
};


template< class Real>
XForm4x4<Real> GetPointStreamScale(vcg::Box3<Real> &bb, float expFact)
{
  qDebug("bbox %f %f %f - %f %f %f ",bb.min[0],bb.min[1],bb.min[2],bb.max[0],bb.max[1],bb.max[2]);
  Real scale = bb.Dim()[bb.MaxDim()] * expFact; 
  Point3m center = bb.Center();
  for( int i=0 ; i<3 ; i++ ) center[i] -= scale/2;
  XForm4x4< Real > tXForm = XForm4x4< Real >::Identity() , sXForm = XForm4x4< Real >::Identity();
  for( int i=0 ; i<3 ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(3,i) = -center[i];
  return sXForm * tXForm;
}


template< class Real >
class MeshModelPointStream : public OrientedPointStreamWithData< Real, Point3m >
{   
  CMeshO &_m;
  size_t _curPos;
public:
  MeshModelPointStream(  CMeshO &m):_m(m),_curPos(0)
  {
    vcg::tri::RequireCompactness(m);
  }
  ~MeshModelPointStream( void ){}
  void reset( void ) { _curPos =0;}
  bool nextPoint( OrientedPoint3D< Real >& pt, Point3m &d)
  {
    if(_curPos>=_m.vn)
      return false;
    Point3m &nn = _m.vert[_curPos].N();
    Point3m tp = _m.Tr * _m.vert[_curPos].P();
    Point4m np = _m.Tr *  Point4m(nn[0],nn[1],nn[2],0);
    
    pt.p[0] = tp[0];
    pt.p[1] = tp[1];
    pt.p[2] = tp[2];
    pt.n[0] = np[0];
    pt.n[1] = np[1];
    pt.n[2] = np[2];
    
    d[0]=Real(_m.vert[_curPos].C()[0]);
    d[1]=Real(_m.vert[_curPos].C()[1]);
    d[2]=Real(_m.vert[_curPos].C()[2]);
    
    ++_curPos;
    return true;
  }

};

template< class Real >
class MeshDocumentPointStream : public OrientedPointStreamWithData< Real, Point3m >
{
  MeshDocument &_md;
  MeshModel *_curMesh;
  size_t _curPos;
  size_t _totalSize;
public:
  MeshDocumentPointStream(  MeshDocument &md):_md(md),_curMesh(0),_curPos(0)
  {
    _totalSize=0;
    MeshModel *m=0;
    do
    {
      m=_md.nextVisibleMesh(m);
      if(m!=0)
      {
        vcg::tri::RequireCompactness(m->cm);
        _totalSize+=m->cm.vn;
      }
    } while(m);
    qDebug("TotalSize %i",_totalSize);
  }
  ~MeshDocumentPointStream( void ){}
  void reset( void ) { _curPos =0; _curMesh=0;}
  bool nextPoint( OrientedPoint3D< Real >& pt, Point3m &d )
  {
    Point3m nn(0,0,0);
//    do
    {
      if((_curMesh==0) || (_curPos >= _curMesh->cm.vn) )
      {
        _curMesh = _md.nextVisibleMesh(_curMesh);
        _curPos = 0;
      }

      if(_curMesh==0)
        return false;
      if(_curPos < _curMesh->cm.vn)
      {
        nn = _curMesh->cm.vert[_curPos].N();
        Point3m tp = _curMesh->cm.Tr * _curMesh->cm.vert[_curPos].P();
        Point4m np = _curMesh->cm.Tr *  Point4m(nn[0],nn[1],nn[2],0);
//        Point3m tp = _curMesh->cm.vert[_curPos].P();
//        Point3m np = nn;
        pt.p[0] = tp[0];
        pt.p[1] = tp[1];
        pt.p[2] = tp[2];
        pt.n[0] = np[0];
        pt.n[1] = np[1];
        pt.n[2] = np[2];
        d[0]=Real(_curMesh->cm.vert[_curPos].C()[0]);
        d[1]=Real(_curMesh->cm.vert[_curPos].C()[1]);
        d[2]=Real(_curMesh->cm.vert[_curPos].C()[2]);
        
        ++_curPos;
      }
    }
    assert(nn!=Point3m(0,0,0));
    return true;
  }

};




template< class Real , int Degree , BoundaryType BType , class Vertex >
int _Execute(OrientedPointStream< Real > *pointStream, Box3m bb, CMeshO &pm, PoissonParam<Real> &pp)
{
	std::cout << " exe start 1" << std::endl;
	/*--------------------------------------------------------------- Start Octree ---------------------------------------------------------------*/
	std::cout << "3 OPENMP threads number : " << pp.ThreadsVal << std::endl;
  typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
  typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
  typedef OrientedPointStreamWithData< Real , Point3D< Real > > PointStreamWithData;
  typedef TransformedOrientedPointStreamWithData< Real , Point3D< Real > > XPointStreamWithData;
  Reset< Real >();
  std::vector< char* > comments;
  
    XForm4x4< Real > xForm = GetPointStreamScale(bb,pp.ScaleVal);
    XForm4x4< Real > iXForm = xForm.inverse();
	DumpOutput2( comments , "Running Screened Poisson Reconstruction (Version 9.0)\n" );
	double startTime = Time();

	OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
	Octree< Real > tree;
	OctreeProfiler< Real > profiler( tree );
	tree.threads = pp.ThreadsVal;
	std::cout <<"4 OPENMP threads number : "<< pp.ThreadsVal << std::endl;
    if( pp.MaxSolveDepthVal<0 ) pp.MaxSolveDepthVal = pp.MaxDepthVal;
    	
	std::cout << " exe start 2" << std::endl;
    
    qDebug("Using %i threads\n",pp.ThreadsVal);
//	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
    if(pp.KernelDepthVal<0) pp.KernelDepthVal =pp.MaxDepthVal-2;
    if( pp.KernelDepthVal>pp.MaxDepthVal )
    {
      printf("kernelDepth cannot be greateer Depth.value\n");
      return false;
    }
    
	int pointCount;

	Real pointWeightSum;
	std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
	std::vector< ProjectiveData< Point3D< Real > , Real > >* sampleData = NULL;
    DensityEstimator* density = NULL;
    SparseNodeData< Point3D< Real > , NORMAL_DEGREE >* normalInfo = NULL;
    Real targetValue = (Real)0.5;
	std::cout << " exe start 3" << std::endl;
	// Read in the samples (and color data)
	{
		profiler.start();
//		PointStream* pointStream;
        
//		char* ext = GetFileExtension( In.value );
//		if( Color.set && Color.value>0 )
//		{
//			sampleData = new std::vector< ProjectiveData< Point3D< Real > , Real > >();
//			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStreamWithData< Real , Point3D< Real > , float , Point3D< unsigned char > >( In.value );
//			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::PlyProperties , 6 , ColorInfo< Real >::ValidPlyProperties );
//			else                                    pointStream = new  ASCIIOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::ReadASCII );
//		}
//		else
//		{
//			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStream< Real , float >( In.value );
//			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStream< Real >( In.value );
//			else                                    pointStream = new  ASCIIOrientedPointStream< Real >( In.value );
//		}
//		delete[] ext;
        sampleData = new std::vector< ProjectiveData< Point3D< Real > , Real > >();        
        XPointStreamWithData _pointStream( xForm , ( PointStreamWithData& )*pointStream );
        pointCount = tree.template init< Point3D< Real > >( _pointStream , pp.MaxDepthVal , pp.ConfidenceFlag , *samples , sampleData );
        
#pragma omp parallel for num_threads( pp.ThreadsVal )
		for( int i=0 ; i<(int)samples->size() ; i++ ) (*samples)[i].sample.data.n *= (Real)-1;

		DumpOutput( "Input Points / Samples: %d / %d\n" , pointCount , samples->size() );
		profiler.dumpOutput2( comments , "# Read input into tree:" );
	}
	std::cout << " exe start 4" << std::endl;
	DenseNodeData< Real , Degree > solution;
	/*--------------------------------------------------------------- Poisson Calcualte ---------------------------------------------------------------*/

	{
		DenseNodeData< Real , Degree > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = pp.MaxSolveDepthVal;

		tree.resetNodeIndices();

		// Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
		{
			profiler.start();
			density = tree.template setDensityEstimator< WEIGHT_DEGREE >( *samples , pp.KernelDepthVal , pp.SamplesPerNodeVal );
			profiler.dumpOutput2( comments , "#   Got kernel density:" );
		}

		// Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
		{
			profiler.start();
			normalInfo = new SparseNodeData< Point3D< Real > , NORMAL_DEGREE >();
			*normalInfo = tree.template setNormalField< NORMAL_DEGREE >( *samples , *density , pointWeightSum , BType==BOUNDARY_NEUMANN );
			profiler.dumpOutput2( comments , "#     Got normal field:" );
		}

		if( !pp.DensityFlag ) delete density , density = NULL;
		std::cout << " exe start 5" << std::endl;
		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			std::vector< int > indexMap;

			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
			tree.template inalizeForBroodedMultigrid< MAX_DEGREE , Degree , BType >( pp.FullDepthVal , typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >( *normalInfo ) , &indexMap );

			if( normalInfo ) normalInfo->remapIndices( indexMap );
			if( density ) density->remapIndices( indexMap );
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}

		// Add the FEM constraints
		{
			profiler.start();
			constraints = tree.template initDenseNodeData< Degree >( );
			tree.template addFEMConstraints< Degree , BType , NORMAL_DEGREE , BType >( FEMVFConstraintFunctor< NORMAL_DEGREE , BType , Degree , BType >( 1. , 0. ) , *normalInfo , constraints , solveDepth );
			profiler.dumpOutput2( comments , "#  Set FEM constraints:" );
		}

		// Free up the normal info [If we don't need it for subseequent iterations.]
		delete normalInfo , normalInfo = NULL;

		// Add the interpolation constraints
		if( pp.PointWeightVal>0 )
		{
			profiler.start();
			iInfo = new InterpolationInfo( tree , *samples , targetValue , pp.AdaptiveExponentVal , (Real)pp.PointWeightVal * pointWeightSum , (Real)0 );
			tree.template addInterpolationConstraints< Degree , BType >( *iInfo , constraints , solveDepth );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
		}

		DumpOutput( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );
		DumpOutput( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename Octree< Real >::SolverInfo solverInfo;
			solverInfo.cgDepth = pp.CGDepthVal , solverInfo.iters = pp.ItersVal , solverInfo.cgAccuracy = pp.CSSolverAccuracyVal , solverInfo.verbose = pp.VerboseFlag , solverInfo.showResidual = pp.ShowResidualFlag , solverInfo.lowResIterMultiplier = std::max< double >( 1. , pp.LowResIterMultiplierVal );
			solution = tree.template solveSystem< Degree , BType >( FEMSystemFunctor< Degree , BType >( 0 , 1. , 0 ) , iInfo , constraints , solveDepth , solverInfo );
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	CoredFileMeshData< Vertex > mesh;

	{
		profiler.start();
		double valueSum = 0 , weightSum = 0;
		typename Octree< Real >::template MultiThreadedEvaluator< Degree , BType > evaluator( &tree , solution , pp.ThreadsVal );
#pragma omp parallel for num_threads( pp.ThreadsVal ) reduction( + : valueSum , weightSum )
		for( int j=0 ; j<samples->size() ; j++ )
		{
			ProjectiveData< OrientedPoint3D< Real > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if( w>0 ) weightSum += w , valueSum += evaluator.value( sample.data.p / sample.weight , omp_get_thread_num() , (*samples)[j].node ) * w;
		}
		Real isoValue = (Real)( valueSum / weightSum );
//		if( samples ) delete samples , samples = NULL;
		profiler.dumpOutput( "Got average:" );
		DumpOutput( "Iso-Value: %e\n" , isoValue );


		/*--------------------------------------------------------------- Mesh reconstruction ---------------------------------------------------------------*/
		profiler.start();
		SparseNodeData< ProjectiveData< Point3D< Real > , Real > , DATA_DEGREE >* colorData = NULL;
		if( sampleData )
		{
			colorData = new SparseNodeData< ProjectiveData< Point3D< Real > , Real > , DATA_DEGREE >();
			*colorData = tree.template setDataField< DATA_DEGREE , false >( *samples , *sampleData , (DensityEstimator*)NULL );
			delete sampleData , sampleData = NULL;
			for( const OctNode< TreeNodeData >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
			{
				ProjectiveData< Point3D< Real > , Real >* clr = (*colorData)( n );
				if( clr ) 
                  (*clr) *= (Real)pow( pp.ColorVal , tree.depth( n ) );
			}
		}
		tree.template getMCIsoSurface< Degree , BType , WEIGHT_DEGREE , DATA_DEGREE >( density , colorData , solution , isoValue , mesh , !pp.LinearFitFlag , !pp.NonManifoldFlag , false /*PolygonMesh.set*/ );
		DumpOutput( "Vertices / Polygons: %d / %d\n" , mesh.outOfCorePointCount()+mesh.inCorePoints.size() , mesh.polygonCount() );
		profiler.dumpOutput2( comments , "#        Got triangles:" );
    }

//        FreePointer( solution );
      
        //cb(90,"Creating Mesh");
        mesh.resetIterator();
        int vm = mesh.outOfCorePointCount()+mesh.inCorePoints.size();
        for(auto pt=mesh.inCorePoints.begin();pt!=mesh.inCorePoints.end();++pt)
        {
          Point3D<Real> pp = iXForm*pt->point;          
           vcg::tri::Allocator<CMeshO>::AddVertex(pm,Point3m(pp[0],pp[1],pp[2]));
           pm.vert.back().Q() = pt->value;
           pm.vert.back().C()[0] = pt->color[0];
           pm.vert.back().C()[1] = pt->color[1];
           pm.vert.back().C()[2] = pt->color[2];
        }
        for (int ii=0; ii < mesh.outOfCorePointCount(); ii++){
          Vertex pt; 
          mesh.nextOutOfCorePoint(pt);
          Point3D<Real> pp = iXForm*pt.point;
          vcg::tri::Allocator<CMeshO>::AddVertex(pm,Point3m(pp[0],pp[1],pp[2]));
          pm.vert.back().Q() = pt.value;
          pm.vert.back().C()[0] = pt.color[0];
          pm.vert.back().C()[1] = pt.color[1];
          pm.vert.back().C()[2] = pt.color[2];
        }
      
        std::vector< CoredVertexIndex > polygon;      
        while(mesh.nextPolygon( polygon ))
        {
          assert(polygon.size()==3);
          int indV[3];
          for( int i=0 ; i<int(polygon.size()) ; i++ )
          {
            if( polygon[i].inCore ) indV[i] = polygon[i].idx;
            else                    indV[i]= polygon[i].idx + int( mesh.inCorePoints.size() );
          }
          vcg::tri::Allocator<CMeshO>::AddFace(pm, &pm.vert[indV[0]], &pm.vert[indV[1]], &pm.vert[indV[2]]);
        }
       // cb(100,"Done");
        
//		if( colorData ) delete colorData , colorData = NULL;


    if( density ) delete density , density = NULL;
	DumpOutput2( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , tree.maxMemoryUsage() );

	return 1;
}



template <class MeshType>
void PoissonClean(MeshType &m, bool scaleNormal, bool cleanFlag)
{
  vcg::tri::UpdateNormal<MeshType>::NormalizePerVertex(m);

  if(cleanFlag) {
    for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if (vcg::SquaredNorm(vi->N()) < std::numeric_limits<float>::min()*10.0)
        vcg::tri::Allocator<MeshType>::DeleteVertex(m,*vi);
 
    for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
        if( fi->V(0)->IsD() || fi->V(1)->IsD() || fi->V(2)->IsD() ) 
          vcg::tri::Allocator<MeshType>::DeleteFace(m,*fi);          
  }
  
  vcg::tri::Allocator<MeshType>::CompactEveryVector(m);
  if(scaleNormal)
  {
    for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
      vi->N() *= vi->Q();
  }
}

bool HasGoodNormal(CMeshO &m)
{
  for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
    if(vcg::SquaredNorm(vi->N()) < std::numeric_limits<float>::min()*10.0) 
      return false; 
  
  return true;
}


//bool importMesh(QString fileName, bool isareload)
//{
//
//
//
//
//	QTime allFileTime;
//	allFileTime.start();
//	//foreach(fileName, fileNameList)
//	{
//		QFileInfo fi(fileName);
//		QString extension = fi.suffix();
//		MeshIOInterface *pCurrentIOPlugin = PM.allKnowInputFormats[extension.toLower()];
//		//pCurrentIOPlugin->setLog(gla->log);
//		if (pCurrentIOPlugin == NULL)
//		{
//			QString errorMsgFormat("Unable to open file:\n\"%1\"\n\nError details: file format " + extension + " not supported.");
//			QMessageBox::critical(this, tr("Meshlab Opening Error"), errorMsgFormat.arg(fileName));
//			return false;
//		}
//
//		RichParameterSet prePar;
//		pCurrentIOPlugin->initPreOpenParameter(extension, fileName, prePar);
//		if (!prePar.isEmpty())
//		{
//			GenericParamDialog preOpenDialog(this, &prePar, tr("Pre-Open Options"));
//			preOpenDialog.setFocus();
//			preOpenDialog.exec();
//		}
//		prePar.join(currentGlobalParams);
//		int mask = 0;
//		//MeshModel *mm= new MeshModel(gla->meshDoc);
//		QFileInfo info(fileName);
//		MeshModel *mm = meshDoc()->addNewMesh(fileName, info.fileName());
//		qb->show();
//		QTime t;
//		t.start();
//		Matrix44m mtr;
//		mtr.SetIdentity();
//		bool open = loadMesh(fileName, pCurrentIOPlugin, mm, mask, &prePar, mtr, isareload);
//		if (open)
//		{
//			GLA()->Logf(0, "Opened mesh %s in %i msec", qUtf8Printable(fileName), t.elapsed());
//			RichParameterSet par;
//			pCurrentIOPlugin->initOpenParameter(extension, *mm, par);
//			if (!par.isEmpty())
//			{
//				GenericParamDialog postOpenDialog(this, &par, tr("Post-Open Processing"));
//				postOpenDialog.setFocus();
//				postOpenDialog.exec();
//				pCurrentIOPlugin->applyOpenParameter(extension, *mm, par);
//			}
//			/*MultiViewer_Container* mv = GLA()->mvc();
//			if (mv != NULL)
//			{
//			for(int glarid = 0;glarid < mv->viewerCounter();++glarid)
//			{
//			GLArea* ar = mv->getViewer(glarid);
//			if (ar != NULL)
//			MLSceneRenderModeAdapter::setupRequestedAttributesAccordingToRenderMode(mm->id(),*ar);
//			}
//			}*/
//		}
//		else
//		{
//			meshDoc()->delMesh(mm);
//			GLA()->Logf(0, "Warning: Mesh %s has not been opened", qUtf8Printable(fileName));
//		}
//	}// end foreach file of the input list
//	GLA()->Logf(0, "All files opened in %i msec", allFileTime.elapsed());
//
//	if (_currviewcontainer != NULL)
//	{
//		_currviewcontainer->resetAllTrackBall();
//		_currviewcontainer->updateAllDecoratorsForAllViewers();
//	}
//	qb->reset();
//	return true;
//}

bool PoissonReconstruction::run(orth::MeshModel& mm,const int depth, const float error_threshold)
{
	std::cout << " program in " << std::endl;
	//vcg::CallBackPos *cb = 0;
	MeshDocument md;
	MeshModel *pm = md.addNewMesh("", "Poisson mesh", false);
	mm.ModelSample(500000);
	orth::MeshModel whole_mesh;

	if (mm.C.size()>0)
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, mm.S.size());
		for (int vertex_index = 0; vertex_index < mm.S.size(); vertex_index++)
		{
			(*vi).P()[0] = (float)(mm.P[mm.S[vertex_index]].x);
			(*vi).P()[1] = (float)(mm.P[mm.S[vertex_index]].y);
			(*vi).P()[2] = (float)(mm.P[mm.S[vertex_index]].z);

			whole_mesh.P.push_back(mm.P[mm.S[vertex_index]]);


			vi->N() = vcg::Point3f((float)(mm.N[mm.S[vertex_index]].x), (float)(mm.N[mm.S[vertex_index]].y), (float)(mm.N[mm.S[vertex_index]].z));

			(*vi).C()[0] = mm.C[mm.S[vertex_index]].x;
			(*vi).C()[1] = mm.C[mm.S[vertex_index]].y;
			(*vi).C()[2] = mm.C[mm.S[vertex_index]].z;
			(*vi).C()[3] = 255;

			++vi;
		}
	}
	else
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, mm.S.size());
		for (int vertex_index = 0; vertex_index < mm.S.size(); vertex_index++)
		{
			(*vi).P()[0] = (float)(mm.P[mm.S[vertex_index]].x);
			(*vi).P()[1] = (float)(mm.P[mm.S[vertex_index]].y);
			(*vi).P()[2] = (float)(mm.P[mm.S[vertex_index]].z);

			whole_mesh.P.push_back(mm.P[mm.S[vertex_index]]);

			vi->N() = vcg::Point3f((float)(mm.N[vertex_index].x), (float)(mm.N[vertex_index].y), (float)(mm.N[vertex_index].z));

			++vi;
		}
	}


	vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(pm->cm);
	std::cout << " changed done " << pm->cm.vert.size()<< std::endl;
	vcg::tri::UpdateBounding<CMeshO>::Box(pm->cm);
	pm->updateDataMask(MeshModel::MM_VERTCOLOR);
	int delVertNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateVertex(pm->cm);
	int delFaceNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateFace(pm->cm);
	vcg::tri::Allocator<CMeshO>::CompactEveryVector(pm->cm);

	
    PoissonParam<Scalarm> pp;
	std::cout << "1 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (depth>0)
	{
		pp.MaxDepthVal = depth;
	}
	else
	{
		pp.MaxDepthVal = 10;
	}
	pp.FullDepthVal = 2;
	pp.CGDepthVal = 0;
	pp.ScaleVal = 1.1;
	pp.SamplesPerNodeVal = 1.5f;
	pp.PointWeightVal = 8.0f;
	pp.ItersVal = 8;
	pp.ConfidenceFlag = false;
	pp.DensityFlag = true;
	pp.CleanFlag = false;
 //   
   bool goodNormal=true, goodColor=true;
	MeshModel *_mm = 0;
	_mm = md.nextVisibleMesh(_mm);
	PoissonClean(_mm->cm, pp.ConfidenceFlag, pp.CleanFlag);
	goodNormal=HasGoodNormal(_mm->cm);
	goodColor = _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
	std::cout << "4 - normal " << goodNormal<<" color - "<<goodColor<< std::endl;

	MeshModel *pm2 = md.addNewMesh("", "Poisson mesh", false);
	pm2->updateDataMask(MeshModel::MM_VERTQUALITY);
	if (goodColor) pm2->updateDataMask(MeshModel::MM_VERTCOLOR);

	MeshModelPointStream<Scalarm> meshStream(_mm->cm);

	std::cout << "2 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (_mm->cm.bbox.max.X()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.X()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(-min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Y()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), min_box_size, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Y()>0)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), 0, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Z()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Z()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), -min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	//_mm->cm.bbox.Set()

	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, _mm->cm.bbox, pm2->cm, pp);

	//{
	//	std::cout << "3" << std::endl;
	//	Box3m bb;
	//	MeshModel *_mm = 0;
	//	int k = 0;
	//	while (_mm = md.nextVisibleMesh(_mm))
	//	{
	//		std::cout << k << std::endl;
	//		bb.Add(_mm->cm.Tr, _mm->cm.bbox);
	//		k++;
	//	}
	//		

	//	MeshDocumentPointStream<Scalarm> documentStream(md);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&documentStream, bb, pm->cm, pp);
	//}
	//else
	//{
	//	std::cout << "4" << std::endl;
	//	MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, md.mm()->cm.bbox, pm->cm, pp);
	//}
	  std::cout << " execute end " << std::endl;
    //pm->UpdateBoxAndNormals();
    //md.setVisible(pm->id(),true);
	//md.setCurrentMesh(pm->id());

	std::cout << " read start " << std::endl;

	int face_index = 0;
	pm2->UpdateBoxAndNormals();

	mm.Clear();
	
	for (auto fi = pm2->cm.face.begin(); fi != pm2->cm.face.end(); ++fi)
	{
		orth::Face face_;
		//fi->cP0
		face_.x = fi->V(0)->Index();
		face_.y = fi->V(1)->Index();
		face_.z = fi->V(2)->Index();


		//for (size_t i = 0; i < 3; i++)
		//{
		//	std::cout << "the " << fi->V(i)->Index() << " point :" << fi->V(i)->P()[0] << ", " << fi->V(i)->P()[1] << ", " << fi->V(i)->P()[2] << "; " << std::endl;
		//}
		//std::cout << " -------------------------------- " << std::endl;

		//std::cout << "the " << fi->Index() << " face :" << face_.x << ", " << face_.y << ", " << face_.z << "; " << std::endl;
		mm.F.push_back(face_);
		face_index++;
	}
	if (goodColor)
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			orth::Color color_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];
			color_.x = pi->C()[0];
			color_.y = pi->C()[1];
			color_.z = pi->C()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
			mm.C.push_back(color_);
		}

	}
	else
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);

		}

	}


	//NearestNeighborSearches nns;
	//vector<unsigned int> query_indexV;
	//vector<double> nearest_distanceV;

	//orth::NearestPointSearch(&whole_mesh, &mm, 8, query_indexV, nearest_distanceV);
	////nns.NearestPointsSearchGPU(&whole_mesh, &mm, 50, query_indexV, nearest_distanceV);

	vector<int> query_indexV;
	vector<float> nearest_distanceV;
	NearestNeighborSearches nns;
	//orth::NearestPointSearch(&mm_temp, &model_to_filt, 8, query_indexV, nearest_distanceV);
	nns.NearestPointsSearchGPU(&whole_mesh, &mm, 50, query_indexV, nearest_distanceV);

	mm.L.resize(mm.P.size(), 0);
	for (int point_index = 0; point_index < mm.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] < error_threshold && nearest_distanceV[point_index] >0.0)
		{
			continue;
		}
		mm.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(mm.P.size(), -1);
	orth::PointColor colors;
	if (goodColor) {

		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				colors.push_back(mm.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < mm.F.size(); face_index++)
	{
		if (mm.L[mm.F[face_index].x] || mm.L[mm.F[face_index].y] || mm.L[mm.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[mm.F[face_index].x];
			f.y = new_point_index[mm.F[face_index].y];
			f.z = new_point_index[mm.F[face_index].z];
			//if (f.x>mm.P.size()|| f.y>mm.P.size()|| f.z>mm.P.size())
			//{
			//	cout << mm.F[face_index].x << "; " << mm.F[face_index].y << "; " << mm.F[face_index].z << endl;
			//	cout << new_point_index[mm.F[face_index].x] << "; " << new_point_index[mm.F[face_index].y] << "; " << new_point_index[mm.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}

	mm.F.swap(faces);
	mm.P.swap(points);
	mm.N.swap(normals);
	if (goodColor)
		mm.C.swap(colors);
	mm.L.swap(labels);

	mm.NormalUpdate();
	mm.S.clear();
	mm.SmallModelFilter(10000);

    return true;

}


bool PoissonReconstruction::run_cpu(orth::MeshModel& mm, const int depth, const float error_threshold)
{
	std::cout << " program in " << std::endl;
	//vcg::CallBackPos *cb = 0;
	MeshDocument md;
	MeshModel *pm = md.addNewMesh("", "Poisson mesh", false);
	mm.ModelSample(500000);
	orth::MeshModel whole_mesh;

	if (mm.C.size()>0)
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, mm.S.size());
		for (int vertex_index = 0; vertex_index < mm.S.size(); vertex_index++)
		{
			(*vi).P()[0] = (float)(mm.P[mm.S[vertex_index]].x);
			(*vi).P()[1] = (float)(mm.P[mm.S[vertex_index]].y);
			(*vi).P()[2] = (float)(mm.P[mm.S[vertex_index]].z);

			whole_mesh.P.push_back(mm.P[mm.S[vertex_index]]);


			vi->N() = vcg::Point3f((float)(mm.N[mm.S[vertex_index]].x), (float)(mm.N[mm.S[vertex_index]].y), (float)(mm.N[mm.S[vertex_index]].z));

			(*vi).C()[0] = mm.C[mm.S[vertex_index]].x;
			(*vi).C()[1] = mm.C[mm.S[vertex_index]].y;
			(*vi).C()[2] = mm.C[mm.S[vertex_index]].z;
			(*vi).C()[3] = 255;

			++vi;
		}
	}
	else
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, mm.S.size());
		for (int vertex_index = 0; vertex_index < mm.S.size(); vertex_index++)
		{
			(*vi).P()[0] = (float)(mm.P[mm.S[vertex_index]].x);
			(*vi).P()[1] = (float)(mm.P[mm.S[vertex_index]].y);
			(*vi).P()[2] = (float)(mm.P[mm.S[vertex_index]].z);

			whole_mesh.P.push_back(mm.P[mm.S[vertex_index]]);

			vi->N() = vcg::Point3f((float)(mm.N[vertex_index].x), (float)(mm.N[vertex_index].y), (float)(mm.N[vertex_index].z));

			++vi;
		}
	}


	vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(pm->cm);
	std::cout << " changed done " << pm->cm.vert.size() << std::endl;
	vcg::tri::UpdateBounding<CMeshO>::Box(pm->cm);
	pm->updateDataMask(MeshModel::MM_VERTCOLOR);
	int delVertNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateVertex(pm->cm);
	int delFaceNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateFace(pm->cm);
	vcg::tri::Allocator<CMeshO>::CompactEveryVector(pm->cm);


	PoissonParam<Scalarm> pp;
	std::cout << "1 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (depth>0)
	{
		pp.MaxDepthVal = depth;
	}
	else
	{
		pp.MaxDepthVal = 10;
	}
	pp.FullDepthVal = 2;
	pp.CGDepthVal = 0;
	pp.ScaleVal = 1.1;
	pp.SamplesPerNodeVal = 1.5f;
	pp.PointWeightVal = 8.0f;
	pp.ItersVal = 8;
	pp.ConfidenceFlag = false;
	pp.DensityFlag = true;
	pp.CleanFlag = false;
	//   
	bool goodNormal = true, goodColor = true;
	MeshModel *_mm = 0;
	_mm = md.nextVisibleMesh(_mm);
	PoissonClean(_mm->cm, pp.ConfidenceFlag, pp.CleanFlag);
	goodNormal = HasGoodNormal(_mm->cm);
	goodColor = _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
	std::cout << "4 - normal " << goodNormal << " color - " << goodColor << std::endl;

	MeshModel *pm2 = md.addNewMesh("", "Poisson mesh", false);
	pm2->updateDataMask(MeshModel::MM_VERTQUALITY);
	if (goodColor) pm2->updateDataMask(MeshModel::MM_VERTCOLOR);

	MeshModelPointStream<Scalarm> meshStream(_mm->cm);

	std::cout << "2 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (_mm->cm.bbox.max.X()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.X()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(-min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Y()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), min_box_size, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Y()>0)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), 0, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Z()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Z()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), -min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	//_mm->cm.bbox.Set()

	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, _mm->cm.bbox, pm2->cm, pp);

	//{
	//	std::cout << "3" << std::endl;
	//	Box3m bb;
	//	MeshModel *_mm = 0;
	//	int k = 0;
	//	while (_mm = md.nextVisibleMesh(_mm))
	//	{
	//		std::cout << k << std::endl;
	//		bb.Add(_mm->cm.Tr, _mm->cm.bbox);
	//		k++;
	//	}
	//		

	//	MeshDocumentPointStream<Scalarm> documentStream(md);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&documentStream, bb, pm->cm, pp);
	//}
	//else
	//{
	//	std::cout << "4" << std::endl;
	//	MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, md.mm()->cm.bbox, pm->cm, pp);
	//}
	std::cout << " execute end " << std::endl;
	//pm->UpdateBoxAndNormals();
	//md.setVisible(pm->id(),true);
	//md.setCurrentMesh(pm->id());

	std::cout << " read start " << std::endl;

	int face_index = 0;
	pm2->UpdateBoxAndNormals();

	mm.Clear();

	for (auto fi = pm2->cm.face.begin(); fi != pm2->cm.face.end(); ++fi)
	{
		orth::Face face_;
		//fi->cP0
		face_.x = fi->V(0)->Index();
		face_.y = fi->V(1)->Index();
		face_.z = fi->V(2)->Index();


		//for (size_t i = 0; i < 3; i++)
		//{
		//	std::cout << "the " << fi->V(i)->Index() << " point :" << fi->V(i)->P()[0] << ", " << fi->V(i)->P()[1] << ", " << fi->V(i)->P()[2] << "; " << std::endl;
		//}
		//std::cout << " -------------------------------- " << std::endl;

		//std::cout << "the " << fi->Index() << " face :" << face_.x << ", " << face_.y << ", " << face_.z << "; " << std::endl;
		mm.F.push_back(face_);
		face_index++;
	}
	if (goodColor)
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			orth::Color color_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];
			color_.x = pi->C()[0];
			color_.y = pi->C()[1];
			color_.z = pi->C()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
			mm.C.push_back(color_);
		}

	}
	else
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);

		}

	}


	NearestNeighborSearches nns;
	vector<unsigned int> query_indexV;
	vector<double> nearest_distanceV;

	orth::NearestPointSearch(&whole_mesh, &mm, 8, query_indexV, nearest_distanceV);
	//nns.NearestPointsSearchGPU(&whole_mesh, &mm, 50, query_indexV, nearest_distanceV);

	mm.L.resize(mm.P.size(), 0);
	for (int point_index = 0; point_index < mm.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] < error_threshold && nearest_distanceV[point_index] >0.0)
		{
			continue;
		}
		mm.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(mm.P.size(), -1);
	orth::PointColor colors;
	if (goodColor) {

		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				colors.push_back(mm.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < mm.F.size(); face_index++)
	{
		if (mm.L[mm.F[face_index].x] || mm.L[mm.F[face_index].y] || mm.L[mm.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[mm.F[face_index].x];
			f.y = new_point_index[mm.F[face_index].y];
			f.z = new_point_index[mm.F[face_index].z];
			//if (f.x>mm.P.size()|| f.y>mm.P.size()|| f.z>mm.P.size())
			//{
			//	cout << mm.F[face_index].x << "; " << mm.F[face_index].y << "; " << mm.F[face_index].z << endl;
			//	cout << new_point_index[mm.F[face_index].x] << "; " << new_point_index[mm.F[face_index].y] << "; " << new_point_index[mm.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}

	mm.F.swap(faces);
	mm.P.swap(points);
	mm.N.swap(normals);
	if (goodColor)
		mm.C.swap(colors);
	mm.L.swap(labels);

	return true;

}


bool PoissonReconstruction::run(vector<orth::MeshModel>& mms, orth::MeshModel& mm, const int depth,const float error_threshold)
{
	std::cout << " program in " << std::endl;
	//vcg::CallBackPos *cb = 0;
	double t,dt;
	t = Time();
	
	MeshDocument md;
	MeshModel *pm = md.addNewMesh("", "Poisson mesh", false);
	orth::MeshModel whole_mesh;


	for (int mm_index = 0; mm_index < mms.size(); mm_index++)
	{
		mms[mm_index].ModelSample(1000000);
	}

	dt = Time();
	printf("									ModelSample --- %9.1f (s)\n", dt - t);
	t = dt;

	unsigned int vertice_number = 0;
	for (int model_index = 0; model_index < mms.size(); model_index++)
	{
		vertice_number += mms[model_index].S.size();
	}
	if (mms[0].C.size()>0)
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, vertice_number);
		for (int model_index = 0; model_index < mms.size(); model_index++)
		{
			for (int vertex_index = 0; vertex_index < mms[model_index].S.size(); vertex_index++)
			{
				(*vi).P()[0] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].x);
				(*vi).P()[1] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].y);
				(*vi).P()[2] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].z);

				whole_mesh.P.push_back(mms[model_index].P[mms[model_index].S[vertex_index]]);

				vi->N() = vcg::Point3f((float)(mms[model_index].N[mms[model_index].S[vertex_index]].x), (float)(mms[model_index].N[mms[model_index].S[vertex_index]].y), (float)(mms[model_index].N[mms[model_index].S[vertex_index]].z));

				(*vi).C()[0] = mms[model_index].C[mms[model_index].S[vertex_index]].x;
				(*vi).C()[1] = mms[model_index].C[mms[model_index].S[vertex_index]].y;
				(*vi).C()[2] = mms[model_index].C[mms[model_index].S[vertex_index]].z;
				(*vi).C()[3] = 255;

				++vi;
			}
		}


	}
	else
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, vertice_number);
		for (int model_index = 0; model_index < mms.size(); model_index++)
		{
			for (int vertex_index = 0; vertex_index < mms[model_index].S.size(); vertex_index++)
			{
				(*vi).P()[0] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].x);
				(*vi).P()[1] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].y);
				(*vi).P()[2] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].z);

				whole_mesh.P.push_back(mms[model_index].P[mms[model_index].S[vertex_index]]);

				vi->N() = vcg::Point3f((float)(mms[model_index].N[vertex_index].x), (float)(mms[model_index].N[vertex_index].y), (float)(mms[model_index].N[vertex_index].z));

				++vi;
			}
		}
	}

	dt = Time();
	printf("									transform data --- %9.1f (s)\n", dt - t);
	t = dt;


	vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(pm->cm);
	std::cout << " changed done " << pm->cm.vert.size() << std::endl;
	vcg::tri::UpdateBounding<CMeshO>::Box(pm->cm);
	pm->updateDataMask(MeshModel::MM_VERTCOLOR);
	int delVertNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateVertex(pm->cm);
	int delFaceNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateFace(pm->cm);
	vcg::tri::Allocator<CMeshO>::CompactEveryVector(pm->cm);


	PoissonParam<Scalarm> pp;
	if (depth>0)
	{
		pp.MaxDepthVal = depth;
	}
	else
	{
		pp.MaxDepthVal = 10;
	}
	pp.FullDepthVal = 2;
	pp.CGDepthVal = 0;
	pp.ScaleVal = 1.1;
	pp.SamplesPerNodeVal = 1.5f;
	pp.PointWeightVal = 8.0f;
	pp.ItersVal = 8;
	pp.ConfidenceFlag = false;
	pp.DensityFlag = true;
	pp.CleanFlag = false;
	//   
	bool goodNormal = true, goodColor = true;
	MeshModel *_mm = 0;
	_mm = md.nextVisibleMesh(_mm);
	PoissonClean(_mm->cm, pp.ConfidenceFlag, pp.CleanFlag);
	goodNormal = HasGoodNormal(_mm->cm);
	goodColor = _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
	std::cout << "4 - normal " << goodNormal << " color - " << goodColor << std::endl;

	MeshModel *pm2 = md.addNewMesh("", "Poisson mesh", false);
	pm2->updateDataMask(MeshModel::MM_VERTQUALITY);
	if (goodColor) pm2->updateDataMask(MeshModel::MM_VERTCOLOR);

	MeshModelPointStream<Scalarm> meshStream(_mm->cm);
	std::cout << "2 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (_mm->cm.bbox.max.X()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.X()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(-min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Y()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), min_box_size, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Y()>0)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), 0, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Z()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Z()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), -min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}

	dt = Time();
	printf("									init reconstruction --- %9.1f (s)\n", dt - t);
	t = dt;

	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, _mm->cm.bbox, pm2->cm, pp);

	dt = Time();
	printf("									reconstruction done --- %9.1f (s)\n", dt - t);
	t = dt;

	//{
	//	std::cout << "3" << std::endl;
	//	Box3m bb;
	//	MeshModel *_mm = 0;
	//	int k = 0;
	//	while (_mm = md.nextVisibleMesh(_mm))
	//	{
	//		std::cout << k << std::endl;
	//		bb.Add(_mm->cm.Tr, _mm->cm.bbox);
	//		k++;
	//	}
	//		

	//	MeshDocumentPointStream<Scalarm> documentStream(md);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&documentStream, bb, pm->cm, pp);
	//}
	//else
	//{
	//	std::cout << "4" << std::endl;
	//	MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, md.mm()->cm.bbox, pm->cm, pp);
	//}
	std::cout << " execute end " << std::endl;
	//pm->UpdateBoxAndNormals();
	//md.setVisible(pm->id(),true);
	//md.setCurrentMesh(pm->id());

	std::cout << " read start " << std::endl;

	int face_index = 0;
	pm2->UpdateBoxAndNormals();

	dt = Time();
	printf("									UpdateBoxAndNormals --- %9.1f (s)\n", dt - t);
	t = dt;
	

	for (auto fi = pm2->cm.face.begin(); fi != pm2->cm.face.end(); ++fi)
	{
		orth::Face face_;
		//fi->cP0
		face_.x = fi->V(0)->Index();
		face_.y = fi->V(1)->Index();
		face_.z = fi->V(2)->Index();


		//for (size_t i = 0; i < 3; i++)
		//{
		//	std::cout << "the " << fi->V(i)->Index() << " point :" << fi->V(i)->P()[0] << ", " << fi->V(i)->P()[1] << ", " << fi->V(i)->P()[2] << "; " << std::endl;
		//}
		//std::cout << " -------------------------------- " << std::endl;

		//std::cout << "the " << fi->Index() << " face :" << face_.x << ", " << face_.y << ", " << face_.z << "; " << std::endl;
		mm.F.push_back(face_);
		face_index++;
	}
	if (goodColor)
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			orth::Color color_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];
			color_.x = pi->C()[0];
			color_.y = pi->C()[1];
			color_.z = pi->C()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
			mm.C.push_back(color_);
		}

	}
	else
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
		}

	}

	std::cout << " recontruction done !!" << std::endl;
	vector<int> query_indexV;
	vector<float> nearest_distanceV;
	NearestNeighborSearches nns;
	//orth::NearestPointSearch(&mm_temp, &model_to_filt, 8, query_indexV, nearest_distanceV);
	nns.NearestPointsSearchGPU(&whole_mesh, &mm, 50, query_indexV, nearest_distanceV);

	dt = Time();
	printf("									NearestPointsSearchGPU --- %9.1f (s)\n", dt - t);
	t = dt;

	mm.L.resize(mm.P.size(), 0);
	for (int point_index = 0; point_index < mm.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] < error_threshold && nearest_distanceV[point_index] >0.0)
		{
			continue;
		}
		mm.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(mm.P.size(), -1);
	orth::PointColor colors;
	if (goodColor) {
		
		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				colors.push_back(mm.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < mm.F.size(); face_index++)
	{
		if (mm.L[mm.F[face_index].x] || mm.L[mm.F[face_index].y] || mm.L[mm.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[mm.F[face_index].x];
			f.y = new_point_index[mm.F[face_index].y];
			f.z = new_point_index[mm.F[face_index].z];
			//if (f.x>mm.P.size()|| f.y>mm.P.size()|| f.z>mm.P.size())
			//{
			//	cout << mm.F[face_index].x << "; " << mm.F[face_index].y << "; " << mm.F[face_index].z << endl;
			//	cout << new_point_index[mm.F[face_index].x] << "; " << new_point_index[mm.F[face_index].y] << "; " << new_point_index[mm.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}


	dt = Time();
	printf("									cut Redundancy data --- %9.1f (s)\n", dt - t);
	t = dt;


	mm.F.swap(faces);
	mm.P.swap(points);
	mm.N.swap(normals);
	if (goodColor)
		mm.C.swap(colors);
	mm.L.swap(labels);



	//for (int mm_index = 0; mm_index < mms.size(); mm_index++)
	//{
	//	mms[mm_index].SmallModelFilter(1000);
	//}
	mm.NormalUpdate();
	mm.S.clear();
	mm.SmallModelFilter(10000);

	dt = Time();
	printf("									SmallModelFilter --- %9.1f (s)\n", dt - t);
	t = dt;

	///*------------------------------------------- model zoom ------------------------------------------*/
	//for (int point_index = 0; point_index < mm.P.size(); point_index++)
	//{
	//	mm.P[point_index] += mm.N[point_index] * zoom_mm_distance*-1;
	//}


	return true;

}

bool PoissonReconstruction::run_cpu(vector<orth::MeshModel>& mms, orth::MeshModel& mm, const int depth, const float error_threshold)
{
	std::cout << " program in " << std::endl;
	//vcg::CallBackPos *cb = 0;

	MeshDocument md;
	MeshModel *pm = md.addNewMesh("", "Poisson mesh", false);
	orth::MeshModel whole_mesh;

	for (int mm_index = 0; mm_index < mms.size(); mm_index++)
	{
		mms[mm_index].SmallModelFilter(1000);
	}

	for (int mm_index = 0; mm_index < mms.size(); mm_index++)
	{
		mms[mm_index].ModelSample(50000);
	}
	unsigned int vertice_number = 0;
	for (int model_index = 0; model_index < mms.size(); model_index++)
	{
		vertice_number += mms[model_index].S.size();
	}
	if (mms[0].C.size()>0)
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, vertice_number);
		for (int model_index = 0; model_index < mms.size(); model_index++)
		{
			for (int vertex_index = 0; vertex_index < mms[model_index].S.size(); vertex_index++)
			{
				(*vi).P()[0] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].x);
				(*vi).P()[1] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].y);
				(*vi).P()[2] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].z);

				whole_mesh.P.push_back(mms[model_index].P[mms[model_index].S[vertex_index]]);

				vi->N() = vcg::Point3f((float)(mms[model_index].N[mms[model_index].S[vertex_index]].x), (float)(mms[model_index].N[mms[model_index].S[vertex_index]].y), (float)(mms[model_index].N[mms[model_index].S[vertex_index]].z));

				(*vi).C()[0] = mms[model_index].C[mms[model_index].S[vertex_index]].x;
				(*vi).C()[1] = mms[model_index].C[mms[model_index].S[vertex_index]].y;
				(*vi).C()[2] = mms[model_index].C[mms[model_index].S[vertex_index]].z;
				(*vi).C()[3] = 255;

				++vi;
			}
		}


	}
	else
	{
		CMeshO::VertexIterator vi = vcg::tri::Allocator<CMeshO>::AddVertices(pm->cm, vertice_number);
		for (int model_index = 0; model_index < mms.size(); model_index++)
		{
			for (int vertex_index = 0; vertex_index < mms[model_index].S.size(); vertex_index++)
			{
				(*vi).P()[0] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].x);
				(*vi).P()[1] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].y);
				(*vi).P()[2] = (float)(mms[model_index].P[mms[model_index].S[vertex_index]].z);

				whole_mesh.P.push_back(mms[model_index].P[mms[model_index].S[vertex_index]]);

				vi->N() = vcg::Point3f((float)(mms[model_index].N[vertex_index].x), (float)(mms[model_index].N[vertex_index].y), (float)(mms[model_index].N[vertex_index].z));

				++vi;
			}
		}
	}


	vcg::tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(pm->cm);
	std::cout << " changed done " << pm->cm.vert.size() << std::endl;
	vcg::tri::UpdateBounding<CMeshO>::Box(pm->cm);
	pm->updateDataMask(MeshModel::MM_VERTCOLOR);
	int delVertNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateVertex(pm->cm);
	int delFaceNum = vcg::tri::Clean<CMeshO>::RemoveDegenerateFace(pm->cm);
	vcg::tri::Allocator<CMeshO>::CompactEveryVector(pm->cm);


	PoissonParam<Scalarm> pp;
	if (depth>0)
	{
		pp.MaxDepthVal = depth;
	}
	else
	{
		pp.MaxDepthVal = 10;
	}
	pp.FullDepthVal = 2;
	pp.CGDepthVal = 0;
	pp.ScaleVal = 1.1;
	pp.SamplesPerNodeVal = 1.5f;
	pp.PointWeightVal = 8.0f;
	pp.ItersVal = 8;
	pp.ConfidenceFlag = false;
	pp.DensityFlag = true;
	pp.CleanFlag = false;
	//   
	bool goodNormal = true, goodColor = true;
	MeshModel *_mm = 0;
	_mm = md.nextVisibleMesh(_mm);
	PoissonClean(_mm->cm, pp.ConfidenceFlag, pp.CleanFlag);
	goodNormal = HasGoodNormal(_mm->cm);
	goodColor = _mm->hasDataMask(MeshModel::MM_VERTCOLOR);
	std::cout << "4 - normal " << goodNormal << " color - " << goodColor << std::endl;

	MeshModel *pm2 = md.addNewMesh("", "Poisson mesh", false);
	pm2->updateDataMask(MeshModel::MM_VERTQUALITY);
	if (goodColor) pm2->updateDataMask(MeshModel::MM_VERTCOLOR);

	MeshModelPointStream<Scalarm> meshStream(_mm->cm);
	std::cout << "2 OPENMP threads number : " << pp.ThreadsVal << std::endl;

	if (_mm->cm.bbox.max.X()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.X()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(-min_box_size, point.Y(), point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Y()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), min_box_size, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Y()>0)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), 0, point.Z());
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.max.Z()<min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}
	if (_mm->cm.bbox.min.Z()>-min_box_size)
	{
		vcg::Point3f point = _mm->cm.bbox.Center();
		Eigen::Vector3f vp(point.X(), point.Y(), -min_box_size);
		point.FromEigenVector(vp);
		_mm->cm.bbox.Add(point);
	}

	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, _mm->cm.bbox, pm2->cm, pp);

	//{
	//	std::cout << "3" << std::endl;
	//	Box3m bb;
	//	MeshModel *_mm = 0;
	//	int k = 0;
	//	while (_mm = md.nextVisibleMesh(_mm))
	//	{
	//		std::cout << k << std::endl;
	//		bb.Add(_mm->cm.Tr, _mm->cm.bbox);
	//		k++;
	//	}
	//		

	//	MeshDocumentPointStream<Scalarm> documentStream(md);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&documentStream, bb, pm->cm, pp);
	//}
	//else
	//{
	//	std::cout << "4" << std::endl;
	//	MeshModelPointStream<Scalarm> meshStream(md.mm()->cm);
	//	_Execute<Scalarm, 2, BOUNDARY_NEUMANN, PlyColorAndValueVertex<Scalarm> >(&meshStream, md.mm()->cm.bbox, pm->cm, pp);
	//}
	std::cout << " execute end " << std::endl;
	//pm->UpdateBoxAndNormals();
	//md.setVisible(pm->id(),true);
	//md.setCurrentMesh(pm->id());

	std::cout << " read start " << std::endl;

	int face_index = 0;
	pm2->UpdateBoxAndNormals();


	for (auto fi = pm2->cm.face.begin(); fi != pm2->cm.face.end(); ++fi)
	{
		orth::Face face_;
		//fi->cP0
		face_.x = fi->V(0)->Index();
		face_.y = fi->V(1)->Index();
		face_.z = fi->V(2)->Index();


		//for (size_t i = 0; i < 3; i++)
		//{
		//	std::cout << "the " << fi->V(i)->Index() << " point :" << fi->V(i)->P()[0] << ", " << fi->V(i)->P()[1] << ", " << fi->V(i)->P()[2] << "; " << std::endl;
		//}
		//std::cout << " -------------------------------- " << std::endl;

		//std::cout << "the " << fi->Index() << " face :" << face_.x << ", " << face_.y << ", " << face_.z << "; " << std::endl;
		mm.F.push_back(face_);
		face_index++;
	}
	if (goodColor)
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			orth::Color color_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];
			color_.x = pi->C()[0];
			color_.y = pi->C()[1];
			color_.z = pi->C()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
			mm.C.push_back(color_);
		}

	}
	else
	{
		for (auto pi = pm2->cm.vert.begin(); pi != pm2->cm.vert.end(); ++pi)
		{
			orth::Point3f point_;
			orth::Normal normal_;
			point_.x = pi->P()[0];
			point_.y = pi->P()[1];
			point_.z = pi->P()[2];
			//std::cout << "the " << pi->Index() << " point :" << point_.x << ", " << point_.y << ", " << point_.z << "; " << std::endl;
			normal_.x = pi->N()[0];
			normal_.y = pi->N()[1];
			normal_.z = pi->N()[2];

			normal_.normalize();
			mm.P.push_back(point_);
			mm.N.push_back(normal_);
		}

	}

	std::cout << " recontruction done !!" << std::endl;
	NearestNeighborSearches nns;

	vector<unsigned int> query_indexV;
	vector<double> nearest_distanceV;

	orth::NearestPointSearch(&whole_mesh, &mm, 8, query_indexV, nearest_distanceV);
	//nns.NearestPointsSearchGPU(&whole_mesh, &mm,50,query_indexV, nearest_distanceV);

	mm.L.resize(mm.P.size(), 0);
	for (int point_index = 0; point_index < mm.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] < error_threshold && nearest_distanceV[point_index] >0.0)
		{
			continue;
		}
		mm.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(mm.P.size(), -1);
	orth::PointColor colors;
	if (goodColor) {

		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				colors.push_back(mm.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < mm.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!mm.L[point_index])
			{
				//cout<< " good ";
				points.push_back(mm.P[point_index]);
				//colors_.push_back(mm.C[point_index]);
				normals.push_back(mm.N[point_index]);
				labels.push_back(mm.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < mm.F.size(); face_index++)
	{
		if (mm.L[mm.F[face_index].x] || mm.L[mm.F[face_index].y] || mm.L[mm.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[mm.F[face_index].x];
			f.y = new_point_index[mm.F[face_index].y];
			f.z = new_point_index[mm.F[face_index].z];
			//if (f.x>mm.P.size()|| f.y>mm.P.size()|| f.z>mm.P.size())
			//{
			//	cout << mm.F[face_index].x << "; " << mm.F[face_index].y << "; " << mm.F[face_index].z << endl;
			//	cout << new_point_index[mm.F[face_index].x] << "; " << new_point_index[mm.F[face_index].y] << "; " << new_point_index[mm.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}

	mm.F.swap(faces);
	mm.P.swap(points);
	mm.N.swap(normals);
	if (goodColor)
		mm.C.swap(colors);
	mm.L.swap(labels);


	///*------------------------------------------- model zoom ------------------------------------------*/
	//for (int point_index = 0; point_index < mm.P.size(); point_index++)
	//{
	//	mm.P[point_index] += mm.N[point_index] * zoom_mm_distance*-1;
	//}

	return true;

}

void PoissonReconstruction::RedundancyFilter(orth::MeshModel &model_target, orth::MeshModel &model_to_filt,double error_threshold)
{

	model_target.S.clear();
	model_target.ModelSample(100000);

	orth::MeshModel mm_temp;
	for (int point_index = 0; point_index < model_target.S.size(); point_index++)
	{
		mm_temp.P.push_back(model_target.P[model_target.S[point_index]]);
	}

	//NearestNeighborSearches nns;
	//vector<int> query_indexV;
	//vector<double> nearest_distanceV;
	//nns.NearestPointsSearchGPU(&mm_temp, &model_to_filt, 50, query_indexV, nearest_distanceV);

	NearestNeighborSearches nns;
	vector<int> query_indexV;
	vector<float> nearest_distanceV;
	//orth::NearestPointSearch(&mm_temp, &model_to_filt, 8, query_indexV, nearest_distanceV);
	nns.NearestPointsSearchGPU(&mm_temp, &model_to_filt, 50, query_indexV, nearest_distanceV);

	model_to_filt.L.clear();
	model_to_filt.L.resize(model_to_filt.P.size(), 0);
	for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] > error_threshold)
		{
			continue;
		}
		model_to_filt.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(model_to_filt.P.size(), -1);
	orth::PointColor colors;
	if (model_to_filt.C.size()== model_to_filt.P.size()) {

		for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!model_to_filt.L[point_index])
			{
				//cout<< " good ";
				points.push_back(model_to_filt.P[point_index]);
				//colors_.push_back(model_to_filt.C[point_index]);
				normals.push_back(model_to_filt.N[point_index]);
				labels.push_back(model_to_filt.L[point_index]);
				colors.push_back(model_to_filt.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!model_to_filt.L[point_index])
			{
				//cout<< " good ";
				points.push_back(model_to_filt.P[point_index]);
				//colors_.push_back(model_to_filt.C[point_index]);
				normals.push_back(model_to_filt.N[point_index]);
				labels.push_back(model_to_filt.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < model_to_filt.F.size(); face_index++)
	{
		if (model_to_filt.L[model_to_filt.F[face_index].x] || model_to_filt.L[model_to_filt.F[face_index].y] || model_to_filt.L[model_to_filt.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[model_to_filt.F[face_index].x];
			f.y = new_point_index[model_to_filt.F[face_index].y];
			f.z = new_point_index[model_to_filt.F[face_index].z];
			//if (f.x>model_to_filt.P.size()|| f.y>model_to_filt.P.size()|| f.z>model_to_filt.P.size())
			//{
			//	cout << model_to_filt.F[face_index].x << "; " << model_to_filt.F[face_index].y << "; " << model_to_filt.F[face_index].z << endl;
			//	cout << new_point_index[model_to_filt.F[face_index].x] << "; " << new_point_index[model_to_filt.F[face_index].y] << "; " << new_point_index[model_to_filt.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}

	model_to_filt.F.swap(faces);
	model_to_filt.P.swap(points);
	model_to_filt.N.swap(normals);
	if (model_to_filt.C.size() == model_to_filt.P.size())
		model_to_filt.C.swap(colors);
	model_to_filt.L.swap(labels);
}

void PoissonReconstruction::RedundancyFilter_cpu(orth::MeshModel &model_target, orth::MeshModel &model_to_filt, double error_threshold)
{
	model_target.S.clear();
	model_target.ModelSample(100000);

	orth::MeshModel mm_temp;
	for (int point_index = 0; point_index < model_target.S.size(); point_index++)
	{
		mm_temp.P.push_back(model_target.P[model_target.S[point_index]]);
	}

	//NearestNeighborSearches nns;
	vector<unsigned int> query_indexV;
	vector<double> nearest_distanceV;
	orth::NearestPointSearch(&mm_temp, &model_to_filt, 8, query_indexV, nearest_distanceV);
	//nns.NearestPointsSearchGPU(&mm_temp, &model_to_filt, 50, query_indexV, nearest_distanceV);

	model_to_filt.L.clear();
	model_to_filt.L.resize(model_to_filt.P.size(), 0);
	for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
	{
		if (nearest_distanceV[point_index] > error_threshold || nearest_distanceV[point_index]<0)
		{
			continue;
		}
		model_to_filt.L[point_index] = 1;
	}

	orth::PointCloudF points;
	orth::Faces faces;
	orth::PointNormal normals;
	orth::PointLabel labels;
	vector<int> new_point_index(model_to_filt.P.size(), -1);
	orth::PointColor colors;
	if (model_to_filt.C.size() == model_to_filt.P.size()) {

		for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!model_to_filt.L[point_index])
			{
				//cout<< " good ";
				points.push_back(model_to_filt.P[point_index]);
				//colors_.push_back(model_to_filt.C[point_index]);
				normals.push_back(model_to_filt.N[point_index]);
				labels.push_back(model_to_filt.L[point_index]);
				colors.push_back(model_to_filt.C[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}
	else
	{
		for (int point_index = 0; point_index < model_to_filt.P.size(); point_index++)
		{
			//cout << "point number "<<point_index;
			if (!model_to_filt.L[point_index])
			{
				//cout<< " good ";
				points.push_back(model_to_filt.P[point_index]);
				//colors_.push_back(model_to_filt.C[point_index]);
				normals.push_back(model_to_filt.N[point_index]);
				labels.push_back(model_to_filt.L[point_index]);
				new_point_index[point_index] = points.size() - 1;

			}
			//cout << endl;
		}
	}


	for (int face_index = 0; face_index < model_to_filt.F.size(); face_index++)
	{
		if (model_to_filt.L[model_to_filt.F[face_index].x] || model_to_filt.L[model_to_filt.F[face_index].y] || model_to_filt.L[model_to_filt.F[face_index].z])
		{
			continue;
		}
		else
		{
			orth::Face f;
			f.x = new_point_index[model_to_filt.F[face_index].x];
			f.y = new_point_index[model_to_filt.F[face_index].y];
			f.z = new_point_index[model_to_filt.F[face_index].z];
			//if (f.x>model_to_filt.P.size()|| f.y>model_to_filt.P.size()|| f.z>model_to_filt.P.size())
			//{
			//	cout << model_to_filt.F[face_index].x << "; " << model_to_filt.F[face_index].y << "; " << model_to_filt.F[face_index].z << endl;
			//	cout << new_point_index[model_to_filt.F[face_index].x] << "; " << new_point_index[model_to_filt.F[face_index].y] << "; " << new_point_index[model_to_filt.F[face_index].z] << endl;
			//	cout << f.x << "; " << f.y << "; " << f.z << endl;
			//}
			faces.push_back(f);
		}
	}

	model_to_filt.F.swap(faces);
	model_to_filt.P.swap(points);
	model_to_filt.N.swap(normals);
	if (model_to_filt.C.size() == model_to_filt.P.size())
		model_to_filt.C.swap(colors);
	model_to_filt.L.swap(labels);
}

bool PoissonReconstruction::MergeWithoutRedundancy(orth::MeshModel &model_target, vector<orth::MeshModel> &models_to_merge, orth::MeshModel &model_merged, double error_threshold)
{


	if (models_to_merge.size() <= 0)
	{
		return false;
	}

	model_merged.Clear();

	for (int group_index = 0; group_index < models_to_merge.size(); group_index++)
	{

		//NearestNeighborSearches nns;
		//vector<int> query_indexV;
		//vector<double> nearest_distanceV;
		//nns.NearestPointsSearchGPU(&models_to_merge[group_index], &model_target, 50, query_indexV, nearest_distanceV);

		NearestNeighborSearches nns;
		vector<int> query_indexV;
		vector<float> nearest_distanceV;
		nns.NearestPointsSearchGPU(&models_to_merge[group_index], &model_target, 50, query_indexV, nearest_distanceV);
		//orth::NearestPointSearch(&models_to_merge[group_index], &model_target, 8, query_indexV, nearest_distanceV);
		//for (int i = 0; i < 1000; i++)
		//{
		//	std::cout << nearest_distanceV[i] << std::endl;
		//}

		model_target.L.clear();
		model_target.L.resize(model_target.P.size(), 0);
		for (int point_index = 0; point_index < model_target.P.size(); point_index++)
		{
			if (nearest_distanceV[point_index] > error_threshold)
			{
				continue;
			}
			model_target.L[point_index] = 1;
		}

		orth::PointCloudF points;
		orth::Faces faces;
		orth::PointNormal normals;
		orth::PointLabel labels;
		vector<int> new_point_index(model_target.P.size(), -1);
		orth::PointColor colors;
		if (model_target.C.size() == model_target.P.size()) {

			for (int point_index = 0; point_index < model_target.P.size(); point_index++)
			{
				//cout << "point number "<<point_index;
				if (!model_target.L[point_index])
				{
					//cout<< " good ";
					points.push_back(model_target.P[point_index]);
					//colors_.push_back(model_target.C[point_index]);
					normals.push_back(model_target.N[point_index]);
					labels.push_back(model_target.L[point_index]);
					colors.push_back(model_target.C[point_index]);
					new_point_index[point_index] = points.size() - 1;

				}
				//cout << endl;
			}
		}
		else
		{
			for (int point_index = 0; point_index < model_target.P.size(); point_index++)
			{
				//cout << "point number "<<point_index;
				if (!model_target.L[point_index])
				{
					//cout<< " good ";
					points.push_back(model_target.P[point_index]);
					//colors_.push_back(model_target.C[point_index]);
					normals.push_back(model_target.N[point_index]);
					labels.push_back(model_target.L[point_index]);
					new_point_index[point_index] = points.size() - 1;

				}
				//cout << endl;
			}
		}


		for (int face_index = 0; face_index < model_target.F.size(); face_index++)
		{
			if (model_target.L[model_target.F[face_index].x] || model_target.L[model_target.F[face_index].y] || model_target.L[model_target.F[face_index].z])
			{
				continue;
			}
			else
			{
				orth::Face f;
				f.x = new_point_index[model_target.F[face_index].x];
				f.y = new_point_index[model_target.F[face_index].y];
				f.z = new_point_index[model_target.F[face_index].z];
				//if (f.x>model_target.P.size()|| f.y>model_target.P.size()|| f.z>model_target.P.size())
				//{
				//	cout << model_target.F[face_index].x << "; " << model_target.F[face_index].y << "; " << model_target.F[face_index].z << endl;
				//	cout << new_point_index[model_target.F[face_index].x] << "; " << new_point_index[model_target.F[face_index].y] << "; " << new_point_index[model_target.F[face_index].z] << endl;
				//	cout << f.x << "; " << f.y << "; " << f.z << endl;
				//}
				faces.push_back(f);
			}
		}

		model_target.F.swap(faces);
		model_target.P.swap(points);
		model_target.N.swap(normals);
		if (model_target.C.size() == model_target.P.size())
			model_target.C.swap(colors);
		model_target.L.swap(labels);

	}

	models_to_merge.push_back(model_target);
	orth::MergeModels(models_to_merge, model_merged);



}

bool PoissonReconstruction::MergeWithoutRedundancy_cpu(orth::MeshModel &model_target, vector<orth::MeshModel> &models_to_merge, orth::MeshModel &model_merged, double error_threshold)
{


	if (models_to_merge.size() <= 0)
	{
		return false;
	}

	model_merged.Clear();

	for (int group_index = 0; group_index < models_to_merge.size(); group_index++)
	{

		//NearestNeighborSearches nns;
		vector<unsigned int> query_indexV;
		vector<double> nearest_distanceV;
		//nns.NearestPointsSearchGPU(&models_to_merge[group_index], &model_target, 50, query_indexV, nearest_distanceV);
		orth::NearestPointSearch(&models_to_merge[group_index], &model_target, 8, query_indexV, nearest_distanceV);
		//for (int i = 0; i < 1000; i++)
		//{
		//	std::cout << nearest_distanceV[i] << std::endl;
		//}

		model_target.L.clear();
		model_target.L.resize(model_target.P.size(), 0);
		for (int point_index = 0; point_index < model_target.P.size(); point_index++)
		{
			if (nearest_distanceV[point_index] > error_threshold|| nearest_distanceV[point_index]<0)
			{
				continue;
			}
			model_target.L[point_index] = 1;
		}

		orth::PointCloudF points;
		orth::Faces faces;
		orth::PointNormal normals;
		orth::PointLabel labels;
		vector<int> new_point_index(model_target.P.size(), -1);
		orth::PointColor colors;
		if (model_target.C.size() == model_target.P.size()) {

			for (int point_index = 0; point_index < model_target.P.size(); point_index++)
			{
				//cout << "point number "<<point_index;
				if (!model_target.L[point_index])
				{
					//cout<< " good ";
					points.push_back(model_target.P[point_index]);
					//colors_.push_back(model_target.C[point_index]);
					normals.push_back(model_target.N[point_index]);
					labels.push_back(model_target.L[point_index]);
					colors.push_back(model_target.C[point_index]);
					new_point_index[point_index] = points.size() - 1;

				}
				//cout << endl;
			}
		}
		else
		{
			for (int point_index = 0; point_index < model_target.P.size(); point_index++)
			{
				//cout << "point number "<<point_index;
				if (!model_target.L[point_index])
				{
					//cout<< " good ";
					points.push_back(model_target.P[point_index]);
					//colors_.push_back(model_target.C[point_index]);
					normals.push_back(model_target.N[point_index]);
					labels.push_back(model_target.L[point_index]);
					new_point_index[point_index] = points.size() - 1;

				}
				//cout << endl;
			}
		}


		for (int face_index = 0; face_index < model_target.F.size(); face_index++)
		{
			if (model_target.L[model_target.F[face_index].x] || model_target.L[model_target.F[face_index].y] || model_target.L[model_target.F[face_index].z])
			{
				continue;
			}
			else
			{
				orth::Face f;
				f.x = new_point_index[model_target.F[face_index].x];
				f.y = new_point_index[model_target.F[face_index].y];
				f.z = new_point_index[model_target.F[face_index].z];
				//if (f.x>model_target.P.size()|| f.y>model_target.P.size()|| f.z>model_target.P.size())
				//{
				//	cout << model_target.F[face_index].x << "; " << model_target.F[face_index].y << "; " << model_target.F[face_index].z << endl;
				//	cout << new_point_index[model_target.F[face_index].x] << "; " << new_point_index[model_target.F[face_index].y] << "; " << new_point_index[model_target.F[face_index].z] << endl;
				//	cout << f.x << "; " << f.y << "; " << f.z << endl;
				//}
				faces.push_back(f);
			}
		}

		model_target.F.swap(faces);
		model_target.P.swap(points);
		model_target.N.swap(normals);
		if (model_target.C.size() == model_target.P.size())
			model_target.C.swap(colors);
		model_target.L.swap(labels);

	}

	models_to_merge.push_back(model_target);
	orth::MergeModels(models_to_merge, model_merged);



}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterScreenedPoissonPlugin)
