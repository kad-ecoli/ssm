// $Id: superpose.cpp,v 1.10 2008/04/22 13:46:09 ccb Exp $
// =================================================================
//
//    09.03.04   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -------------------------------------------------------------------
//
//  **** Module  :  Superposer <implementation>
//       ~~~~~~~~~
//  **** Project :  Structure alignment in 3D
//       ~~~~~~~~~
//
//  E. Krissinel, 2002-2004
//
// =================================================================
//

#ifndef  __STRING_H
#include <string.h>
#endif

#ifndef  __MATH_H
#include <math.h>
#endif

#ifndef  __SSM_Align__
#include "ssm_align.h"
#endif


void printInstructions ( pstr argv0 )  {

  printf (
  "\n"
  " Secondary Structure Superposition\n"
  " ---------------------------------\n"
  "\n"
  " USAGE:\n"
  "\n"
  "%s s1.pdb [-s CID1] s2.pdb [-s CID2] [ foo_out.pdb ]\n"
  "\n"
  "where  s1.pdb is the Query structure to which transformation applied,\n"
  "       s2.pdb is the fixed Target structure,\n"
  "       [-s CID1/2] are optional selection strings in MMDB convention, and\n"
  "       [foo_out.pdb] is optional output file.\n",argv0
   );

}


void printccp4rot ( mat44 m )
{
  double w_, x_, y_, z_;
  double d = 180.0/3.14159265359;
  double tr = m[0][0] + m[1][1] + m[2][2] + 1.0;
  // check the diagonal
  if ( tr > 1.0e-8 ) {
    double s( sqrt(tr) );
    w_ = s * 0.5;
    s = 0.5 / s;
    x_ = s * ( m[2][1] - m[1][2] );
    y_ = s * ( m[0][2] - m[2][0] );
    z_ = s * ( m[1][0] - m[0][1] );
  } else {
    if ( m[0][0] > m[1][1] && m[0][0] > m[2][2] ) {
      double s( sqrt(1.0 + m[0][0] - m[1][1] - m[2][2] ) );
      x_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m[2][1] - m[1][2] );
      y_ = s * ( m[0][1] + m[1][0] );
      z_ = s * ( m[0][2] + m[2][0] );
    } else if ( m[1][1] > m[2][2] ) {
      double s( sqrt(1.0 + m[1][1] - m[2][2] - m[0][0] ) );
      y_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m[0][2] - m[2][0] );
      z_ = s * ( m[1][2] + m[2][1] );
      x_ = s * ( m[1][0] + m[0][1] );
    } else {
      double s( sqrt(1.0 + m[2][2] - m[0][0] - m[1][1] ) );
      z_ = 0.5 * s;
      if ( s != 0.0 ) s = 0.5 / s;
      w_ = s * ( m[1][0] - m[0][1] );
      x_ = s * ( m[2][0] + m[0][2] );
      y_ = s * ( m[2][1] + m[1][2] );
    }
  }
  double om, ph, ka, al, be, ga;
  om = ph = ka = 0.0;
  if ( fabs(w_) < 0.999999 ) {
    double r = sqrt( x_*x_ + y_*y_ );
    om = d*atan2( r, z_ );
    if ( r > 0.000001 ) ph = d*atan2( y_, x_ );
    ka = d*2.0*acos( w_ );
  }
  double ca, cb, cg, sa, sb, sg;
  cb = 1.0 - 2.0 * (x_*x_ + y_*y_);
  sb = 2.0 * sqrt( (x_*x_ + y_*y_) * (w_*w_ + z_*z_) );
  if ( sb > 0.0001 ) {
    ca = 2.0 * (x_*z_ + w_*y_);
    sa = 2.0 * (y_*z_ - w_*x_);
    cg = 2.0 * (w_*y_ - x_*z_);
    sg = 2.0 * (y_*z_ + w_*x_);
  } else {
    ca = 1.0;
    sa = 0.0;
    cg = cb;
    sg = 2.0*(y_*z_ + w_*x_);
  }
  al = d*atan2(sa,ca);
  be = d*atan2(sb,cb);
  ga = d*atan2(sg,cg);

  printf( "\nCCP4 format rotation-translation operator\n" );
  printf( "Polar angles (omega,phi,kappa) : %9.3f %9.3f %9.3f\n", om, ph, ka );
  printf( "Euler angles (alpha,beta,gamma): %9.3f %9.3f %9.3f\n", al, be, ga );
  printf( "Orthogonal translation (/Angst): %9.3f %9.3f %9.3f\n", m[0][3],m[1][3],m[2][3] );
}



int  readCoorFile ( pstr FName, RPCMMDBManager MMDB )  {
char S[500];
int  rc,lcount;

  if (!MMDB)  MMDB = new CMMDBManager();

  MMDB->SetFlag ( MMDBF_PrintCIFWarnings       |
                  MMDBF_IgnoreNonCoorPDBErrors |
                  MMDBF_IgnoreDuplSeqNum );

  rc = MMDB->ReadCoorFile ( FName );

  if (rc) {
    printf ( " ***** ERROR #%i READ:\n\n %s\n\n",rc,GetErrorDescription(rc) );
    MMDB->GetInputBuffer ( S,lcount );
    if (lcount>=0) 
      printf ( "       LINE #%i:\n%s\n\n",lcount,S );
    else if (lcount==-1)
      printf ( "       CIF ITEM: %s\n\n",S );
    delete MMDB;
    MMDB = NULL;
    return 1;
  } else  {
    switch (MMDB->GetFileType())  {
      case MMDB_FILE_PDB    : printf ( " PDB"         );  break;
      case MMDB_FILE_CIF    : printf ( " mmCIF"       );  break;
      case MMDB_FILE_Binary : printf ( " MMDB binary" );  break;
      default : printf ( " Unknown (report as a bug!)" );
    }
    printf ( " file %s has been read in.\n",FName );
  }

  return 0;

}

int selectAtoms ( PCMMDBManager M, char ** argv, int & argNo,
                  int & selHnd )  {
int nSel;
  selHnd = 0;
  if (!strcasecmp(argv[argNo],"-s"))  {
    argNo++;
    selHnd = M->NewSelection();
    M->Select ( selHnd,STYPE_ATOM,argv[argNo],SKEY_NEW );
    nSel = M->GetSelLength ( selHnd );
    if (nSel<=0)  {
      printf ( " *** Selection string '%s' does not cover any atoms.\n",
               argv[argNo] );
      return 1;
    }
    printf ( " ... %i atoms selected using CID '%s'\n",nSel,argv[argNo] );
    argNo++;
  } else {
    selHnd = M->NewSelection();
    M->Select ( selHnd,STYPE_ATOM,"*",SKEY_NEW );
    nSel = M->GetSelLength ( selHnd );
    if (nSel<=0)  {
      printf ( " *** No atoms in PDB file!\n" );
      return 1;
    }
    printf ( " ... %i atoms read from PDB file\n",nSel );
  }    
  return 0;
}

int main ( int argc, char ** argv, char ** env )  {
CFile         f;
PCMMDBManager M1,M2;
PCSSMAlign    SSMAlign;
pstr          name1,name2,fileout;
int           argNo,selHnd1,selHnd2,rc,n;

  if (argc<=1)  {
    printInstructions ( argv[0] );
    return 1;
  }
  if ((!strcmp(argv[1],"-?")) || (!strcasecmp(argv[1],"-help")) ||
      (!strcasecmp(argv[1],"--help")))  {
    printInstructions ( argv[0] );
    return 1;
  }

  InitMatType();
  InitSSGraph();


  M1    = NULL;
  name1 = NULL;
  M2    = NULL;
  name2 = NULL;
  fileout = NULL;
  argNo = 1;

  CreateCopy ( name1,argv[argNo] );
  if (readCoorFile(argv[argNo++],M1))
    return 2;
  // If argument at argNo is "-s", the selection string is read
  // and argNo updated. Else argNo is unchanged.
  if (selectAtoms(M1,argv,argNo,selHnd1))  {
    delete M1;
    return 3;
  }

  CreateCopy ( name2,argv[argNo] );
  if (readCoorFile(argv[argNo++],M2))  {
    delete M1;
    return 4;
  }
  // n = 0 means no arguments left, but we still want to call
  // selectAtoms to set default selection string
  n = 0;
  if (argNo<argc) n = argNo;
  // If argument at argNo is "-s", the selection string is read
  // and n updated. Else n is unchanged.
  if (selectAtoms(M2,argv,n,selHnd2))  {
    delete M2;
    return 5;
  }
  if (argNo<argc) argNo = n;

  // get optional output file
  if (argNo<argc) {
    CreateCopy ( fileout,argv[argNo] );
    printf ( "\n Transformed coordinates will be written to file %s \n",fileout );
  }

  SSMAlign = new CSSMAlign();
  rc = SSMAlign->Align ( M1,M2,SSMP_Normal,CSSC_Flexible,selHnd1,selHnd2 );

  if (rc)  {
    switch (rc)  {
      case SSM_noHits :
         printf ( " *** secondary structure does not match.\n" );
         break;
      case SSM_noSPSN :
         printf ( " *** structures are too remote.\n" );
         break;
      case SSM_noGraph :
         printf ( " *** can't make graph for %s.\n",name1 );
         break;
      case SSM_noVertices :
         printf ( " *** empty graph for %s.\n",name1 );
         break;
      case SSM_noGraph2 :
         printf ( " *** can't make graph for %s.\n",name2 );
         break;
      case SSM_noVertices2 :
         printf ( " *** empty graph for %s.\n",name2 );
         break;
      default :
        printf ( " *** undocumented return code %i.\n",rc );
    }
  } else  {
    // MrBUMP now uses RMSD line so be careful with formatting
    printf ( "\n"
             " Query      %s\n"
             " and Target %s\n"
             " have been superposed. Transformation matrix (to be applied\n"
             " to %s) is\n\n"
             "        Rx         Ry         Rz           T\n"
             " %10.4f %10.4f %10.4f   %10.4f\n"
             " %10.4f %10.4f %10.4f   %10.4f\n"
             " %10.4f %10.4f %10.4f   %10.4f\n\n"
             " at RMSD = %10.3f and alignment length %i\n",
             name1,name2,name1,
             SSMAlign->TMatrix[0][0],SSMAlign->TMatrix[0][1],
             SSMAlign->TMatrix[0][2],SSMAlign->TMatrix[0][3],
             SSMAlign->TMatrix[1][0],SSMAlign->TMatrix[1][1],
             SSMAlign->TMatrix[1][2],SSMAlign->TMatrix[1][3],
             SSMAlign->TMatrix[2][0],SSMAlign->TMatrix[2][1],
             SSMAlign->TMatrix[2][2],SSMAlign->TMatrix[2][3],
             SSMAlign->rmsd,SSMAlign->nalgn );

    f.assign ( "stdout" );
    f.rewrite();
    PrintSSMAlignTable ( f,M1,M2,SSMAlign );

    f.shut();

    // if output file requested, apply transform to first
    // input file (NB all file, not just requested selection)
    // and output
    if ( fileout ) {
      M1->ApplyTransform ( SSMAlign->TMatrix );
      M1->WritePDBASCII ( fileout );
    }

    printccp4rot( SSMAlign->TMatrix );
  }

  if (M1)       delete M1;
  if (name1)    delete name1;
  if (M2)       delete M2;
  if (name2)    delete name2;
  if (SSMAlign) delete SSMAlign;

  return 0;

}

