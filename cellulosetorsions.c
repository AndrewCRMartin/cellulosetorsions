/************************************************************************/
/**

   \file       cellulosetorsions.c
   
   \version    V1.1
   \date       29.11.17
   \brief      Calculate torsion angles for a PDB file of cellulose 
               conformations
   
   \copyright  (c) Dr. Andrew C. R. Martin 2017
   \author     Dr. Andrew C. R. Martin
   \par
               Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Calculate torsion angles from a PDB file

**************************************************************************

   Usage:
   ======
   NOTE! If the executable is called 'torsions' rather than 'cellulosetorsions'
   this has the effect of setting the '-o' (old style) flag by default.
   Consequently a symbolic link to the program can be used to obtain
   old-style output for backwards compatibility

**************************************************************************

   Revision History:
   =================

-  V1.0  28.11.17 Original
-  V1.1  29.11.17 Does omega as well and creates distribution. Also allows
                  a set of atoms to be specified

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF     512
#define ERROR_VALUE 9999.0
#define MAXDIST     361

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *Radians, BOOL *doDistrib, BOOL *doOther,
                  char *atom1, char *atom2, char *atom3, char *atom4); 
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians);
void CalculateAndDisplayTorsions(FILE *out, PDB *fullpdb, 
                                 BOOL Radians, BOOL doDistrib);
void CalculateAndDisplayOtherTorsion(FILE *out, PDB *pdb, BOOL Radians, 
                                     BOOL doDistrib, char *atom1, 
                                     char *atom2,char *atom3, 
                                     char *atom4);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for converting a PDB file to torsions.

-  28.11.17 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    inFile[MAXBUFF],
           outFile[MAXBUFF];
   int     natoms;
   PDB     *pdb;
   BOOL    Radians    = FALSE,
           doDistrib  = FALSE;

   BOOL doOther = FALSE;
   
   char atom1[8], atom2[8], atom3[8], atom4[8]; 

   if(ParseCmdLine(argc, argv, inFile, outFile, &Radians, &doDistrib,
                   &doOther, atom1, atom2, atom3, atom4))
   {
      if(blOpenStdFiles(inFile, outFile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))!=NULL)
         {
            if(doOther)
            {
               CalculateAndDisplayOtherTorsion(out, pdb, Radians, doDistrib,
                                               atom1, atom2, atom3, atom4);
            }
            else
            {
               CalculateAndDisplayTorsions(out, pdb, Radians, doDistrib);
            }
         }
         else
         {
            fprintf(stderr,"cellulosetorsions: Error - no atoms read \
rom PDB file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"cellulosetorsions: Error - unable to open \
input or output file\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
void CalculateAndDisplayTorsions(FILE *out, PDB *fullpdb, BOOL Radians,
                                 BOOL doDistrib)
{
   PDB *res     = NULL,
       *nextRes = NULL;
   int phiDist[MAXDIST], psiDist[MAXDIST], omegaDist[MAXDIST], i;

   for(i=0; i<MAXDIST; i++)
      phiDist[i] = psiDist[i] = omegaDist[i] = 0;
   
   res = fullpdb;

   if(doDistrib)
   {
      fprintf(out, "#        PHI Count    PSI Count  OMEGA Count\n");
      fprintf(out, "#Angle O5-C1-O4-C4  C1-O4-C4-C5  O5-C5-C6-O6\n");
   }
   else
   {
      fprintf(out, "#        PHI          PSI        OMEGA\n");
      fprintf(out, "#O5-C1-O4-C4  C1-O4-C4-C5  O5-C5-C6-O6\n");
   }
   
   while(res!=NULL)
   {
      BOOL OK;
      int  i;
      PDB *p, 
          *phiAtoms[4],
          *psiAtoms[4],
          *omegaAtoms[4];

      nextRes = blFindNextResidue(res);

      /* Blank the arrays                                               */
      for(i=0; i<4; i++)
         phiAtoms[i] = psiAtoms[i] = omegaAtoms[i] = NULL;
      
      /* Find the atoms for phi and psi in this residue                 */
      for(p=res; p!=nextRes; NEXT(p))
      {
         if(!strncmp(p->atnam, "O5  ", 4))
            phiAtoms[0] = p;
         else if(!strncmp(p->atnam, "C1  ", 4))
            phiAtoms[1] = p;
         else if(!strncmp(p->atnam, "O4  ", 4))
            phiAtoms[2] = p;
         else if(!strncmp(p->atnam, "C4  ", 4))
            phiAtoms[3] = p;
         
         if(!strncmp(p->atnam, "C1  ", 4))
            psiAtoms[0] = p;
         else if(!strncmp(p->atnam, "O4  ", 4))
            psiAtoms[1] = p;
         else if(!strncmp(p->atnam, "C4  ", 4))
            psiAtoms[2] = p;
         else if(!strncmp(p->atnam, "C5  ", 4))
            psiAtoms[3] = p;

         if(!strncmp(p->atnam, "O5  ", 4))
            omegaAtoms[0] = p;
         else if(!strncmp(p->atnam, "C5  ", 4))
            omegaAtoms[1] = p;
         else if(!strncmp(p->atnam, "C6  ", 4))
            omegaAtoms[2] = p;
         else if(!strncmp(p->atnam, "O6  ", 4))
            omegaAtoms[3] = p;
      }
      
      OK = TRUE;
      for(i=0; i<4; i++)
      {
         if((phiAtoms[i]   == NULL) || 
            (psiAtoms[i]   == NULL) || 
            (omegaAtoms[i] == NULL))
         {
            OK = FALSE;
            break;
         }
      }
      
      if(!OK)
      {
         fprintf(stderr, "Atom missing in residue %s.%d%s\n",
                 res->chain, res->resnum, res->insert);
      }
      else
      {
         REAL phi, psi, omega;
         phi   = CalcTorsion(phiAtoms[0],   phiAtoms[1],
                             phiAtoms[2],   phiAtoms[3], Radians);
         psi   = CalcTorsion(psiAtoms[0],   psiAtoms[1],
                             psiAtoms[2],   psiAtoms[3], Radians);
         omega = CalcTorsion(omegaAtoms[0], omegaAtoms[1],
                             omegaAtoms[2], omegaAtoms[3], Radians);
         if(doDistrib)
         {
            phiDist[180+(int)(phi+0.5)]++;
            psiDist[180+(int)(psi+0.5)]++;
            omegaDist[180+(int)(omega+0.5)]++;
         }
         else
         {
            fprintf(out, " %11.2f  %11.2f  %11.2f\n", phi, psi, omega);
         }
      }
      
      /* Step to next residue                                           */
      res=nextRes;
   }

   if(doDistrib)
   {
      for(i=0; i<MAXDIST; i++)
      {
         fprintf(out, "%6d %11d %11d %11d\n", i-180, phiDist[i], 
                 psiDist[i], omegaDist[i]);
      }
   }
}


/************************************************************************/
/*>REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
   ------------------------------------------------------------------
*//**
   \param[in]    *p1      Pointer to PDB record
   \param[in]    *p2      Pointer to PDB record
   \param[in]    *p3      Pointer to PDB record
   \param[in]    *p4      Pointer to PDB record
   \param[in]    Radians  Output radians rather than degrees?
   \return                Torsion angle

   Wraps around the blPhi() function to take (and check) PDB pointers
   rather than coordinates, handle conversion from radians, etc.

- 27.11.14 Original   By: ACRM
*/
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
{
   REAL tor;
   
   if((p1==NULL)||(p2==NULL)||(p3==NULL)||(p4==NULL))
   {
      return(ERROR_VALUE);
   }

   tor = blPhi(p1->x, p1->y, p1->z,
               p2->x, p2->y, p2->z,
               p3->x, p3->y, p3->z,
               p4->x, p4->y, p4->z);

   if(!Radians)
      tor *= 180/PI;

   return(tor);
}

/************************************************************************/
void CalculateAndDisplayOtherTorsion(FILE *out, PDB *fullpdb, 
                                     BOOL Radians, BOOL doDistrib,
                                     char *atom1, char *atom2,
                                     char *atom3, char *atom4)
{
   PDB *res     = NULL,
       *nextRes = NULL;
   int angleDist[MAXDIST], i;

   for(i=0; i<MAXDIST; i++)
      angleDist[i] = 0;
   
   res = fullpdb;

   if(doDistrib)
   {
      fprintf(out, "#      Count\n");
      fprintf(out, "#Angle %s-%s-%s-%s\n", atom1, atom2, atom3, atom4);
   }
   else
   {
      fprintf(out, "#      Angle\n");
      fprintf(out, "#%s-%s-%s-%s\n", atom1, atom2, atom3, atom4);
   }

   PADMINTERM(atom1, 4);
   PADMINTERM(atom2, 4);
   PADMINTERM(atom3, 4);
   PADMINTERM(atom4, 4);
   
   while(res!=NULL)
   {
      BOOL OK;
      int  i;
      PDB *p, 
          *atoms[4];

      nextRes = blFindNextResidue(res);

      /* Blank the arrays                                               */
      for(i=0; i<4; i++)
         atoms[i] = NULL;
      
      /* Find the atoms for phi and psi in this residue                 */
      for(p=res; p!=nextRes; NEXT(p))
      {
         if(!strncmp(p->atnam, atom1, 4))
            atoms[0] = p;
         else if(!strncmp(p->atnam, atom2, 4))
            atoms[1] = p;
         else if(!strncmp(p->atnam, atom3, 4))
            atoms[2] = p;
         else if(!strncmp(p->atnam, atom4, 4))
            atoms[3] = p;
      }
      
      OK = TRUE;
      for(i=0; i<4; i++)
      {
         if(atoms[i] == NULL)
         {
            OK = FALSE;
            break;
         }
      }
      
      if(!OK)
      {
         fprintf(stderr, "Atom missing in residue %s.%d%s\n",
                 res->chain, res->resnum, res->insert);
      }
      else
      {
         REAL angle;
         angle = CalcTorsion(atoms[0],   atoms[1],
                             atoms[2],   atoms[3], Radians);
         if(doDistrib)
         {
            angleDist[180+(int)(angle+0.5)]++;
         }
         else
         {
            fprintf(out, " %11.2f\n", angle);
         }
      }
      
      /* Step to next residue                                           */
      res=nextRes;
   }

   if(doDistrib)
   {
      for(i=0; i<MAXDIST; i++)
      {
         fprintf(out, "%6d %11d\n", i-180, angleDist[i]);
      }
   }
   
}



/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *Radians, BOOL *doDistrib, BOOL *doOther,
                     char *atom1, char *atom2, char *atom3, char *atom4)
   ---------------------------------------------------------------------
*//**

   \param[in]     argc         Argument count
   \param[in]     **argv       Argument array
   \param[out]    *infile      Input file (or blank string)
   \param[out]    *outfile     Output file (or blank string)
   \param[out]    *Radians     Output radians rather than degrees
   \param[out]    *doDistrib   Show the distribution instead of angles
   \param[out]    *doOther     Do a specified set of atoms
   \param[out]    *atom1       Atom1
   \param[out]    *atom2       Atom2
   \param[out]    *atom3       Atom3
   \param[out]    *atom4       Atom4

   \return                     Success?

   Parse the command line
   
-  28.11.17 V1.0
-  29.11.17 Added -d/doDistrib, -a/doOther/atom1/atom2/atom3/atom4
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *Radians, BOOL *doDistrib, BOOL *doOther,
                  char *atom1, char *atom2, char *atom3, char *atom4)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'r':
            *Radians = TRUE;
            break;
         case 'd':
            *doDistrib = TRUE;
            break;
         case 'a':
            *doOther = TRUE;
            argc--; argv++;
            if(argc < 4)
               return(FALSE);
            strncpy(atom1, argv[0], 8);
            argc--; argv++;
            strncpy(atom2, argv[0], 8);
            argc--; argv++;
            strncpy(atom3, argv[0], 8);
            argc--; argv++;
            strncpy(atom4, argv[0], 8);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are <= 2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         /* Copy the second to outfile                                  */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Displays a usage message

-  28.11.17 Original   By: ACRM
-  29.11.17 V1.1 Added -d and -a
*/
void Usage(void)
{
   fprintf(stderr,"\ncellulosetorsions V1.1 (c) 2017 Andrew Martin, \
UCL.\n");
   fprintf(stderr,"\nUsage: cellulosetorsions [-h][-r][-d][-a atom1 \
atom2 atom3 atom4]\n");
   fprintf(stderr,"                         [in.pdb [out.tor]]\n");
   fprintf(stderr,"       -h   This help message\n");
   fprintf(stderr,"       -r   Give results in radians\n");
   fprintf(stderr,"       -d   Show the distribution in 1 degree bins\n");
   fprintf(stderr,"       -a   Instead of Phi,Psi,Omega use the \
specified atoms\n");

   fprintf(stderr,"\nGenerates a set of phi (O5-C1-O4-C4), psi \
(C1-O4-C4-C5) and omega\n");
   fprintf(stderr,"(O5-C5-C6-O6) torsion angles from a PDB file \
containing conformations\n");
   fprintf(stderr,"of cellulose.\n");

   fprintf(stderr,"\nI/O is through stdin/stdout if unspecified.\n\n");
}


