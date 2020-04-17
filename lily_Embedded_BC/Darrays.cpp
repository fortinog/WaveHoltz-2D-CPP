#include <iostream>

#include "Darrays.h"

using namespace std;

// Default value of array ordering
bool Darray1::m_corder = false;
bool Darray2::m_corder = false;

// Constructors
//-----------------------------------------------------------------------
Darray1::Darray1( int nc, int ibeg, int iend)
{
    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
Darray1::Darray1(int ibeg, int iend)
{
    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
Darray1::Darray1()
{
    m_nc = m_ib = m_ie = 0;
    m_data = NULL;
}

//
void Darray1::define( int nc, int ibeg, int iend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}

void Darray1::define(int ibeg, int iend)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
        m_data = new double[m_nc*m_ni];
    else
        m_data = NULL;
    define_offsets();
}



//-----------------------------------------------------------------------
void Darray1::set_value(double scalar)
{
    for( size_t i = 0 ; i < m_npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Darray1::define_offsets()
{
    m_npts = static_cast<size_t>(m_ni)*m_nc;
    if( m_corder )
    {
            // (i,c)=i-ib+ni*(c-1)
        m_base = -m_ib-m_ni;
        m_offc = m_ni;
        m_offi = 1;
    }
    else
    {
            // (c,i)=c-1+nc*(i-ib)
        m_base = -1-m_nc*m_ib;
        m_offc = 1;
        m_offi = m_nc;
    }
}

//-----------------------------------------------------------------------
void Darray1::copy( const Darray1& u )
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = u.m_nc;
    m_ib = u.m_ib;
    m_ie = u.m_ie;
    m_ni = m_ie-m_ib+1;
    if( m_nc*m_ni > 0 )
    {
        m_data = new double[m_nc*m_ni];
        for( int i = 0 ; i < m_nc*m_ni ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}


//-------------------------
// 2D
//-------------------------

// Constructors
//-----------------------------------------------------------------------
Darray2::Darray2( int nc, int ibeg, int iend, int jbeg, int jend)
{
    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    m_jb = jbeg;
    m_je = jend;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
Darray2::Darray2(int ibeg, int iend, int jbeg, int jend)
{
    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    m_jb = jbeg;
    m_je = jend;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
Darray2::Darray2()
{
    m_nc = m_ib = m_ie = m_jb = m_je = 0;
    m_data = NULL;
}

//
void Darray2::define( int nc, int ibeg, int iend, int jbeg, int jend)
{
    if( m_data != NULL )
        delete[] m_data;
    
    m_nc = nc;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    m_jb = jbeg;
    m_je = jend;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}

void Darray2::define(int ibeg, int iend, int jbeg, int jend)
{
    if( m_data != NULL )
        delete[] m_data;
    
    m_nc = 1;
    m_ib = ibeg;
    m_ie = iend;
    m_ni = m_ie-m_ib+1;
    m_jb = jbeg;
    m_je = jend;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_ni*m_nj > 0 )
        m_data = new double[m_nc*m_ni*m_nj];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------
void Darray2::set_value(double scalar)
{
    for( size_t i = 0 ; i < m_npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Darray2::define_offsets()
{
    m_npts = static_cast<size_t>(m_ni)*m_nj*m_nc;
    if( m_corder )
    {
      // (i,j,k,c)=i-ib+ni*(j-jb)+ni*nj*(c-1)
        m_base = -m_ib-m_ni*m_jb;
        m_offc = m_ni*m_nj;
        m_offj = m_ni;
        m_offi = 1;
    }
    else
    {
            // (c,i,j,k)=c-1+nc*(i-ib)+nc*ni*(j-jb);
        m_base = -1-m_nc*m_ib-m_nc*m_ni*m_jb;
        m_offc = 1;
        m_offi = m_nc;
        m_offj = m_nc*m_ni;
    }
}

//-----------------------------------------------------------------------
void Darray2::copy( const Darray2& u )
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = u.m_nc;
    m_ib = u.m_ib;
    m_ie = u.m_ie;
    m_ni = m_ie-m_ib+1;
    m_jb = u.m_jb;
    m_je = u.m_je;
    m_nj = m_je-m_jb+1;
    if( m_nc*m_nj*m_ni > 0 )
    {
        m_data = new double[m_nc*m_ni*m_nj];
        for( int i = 0 ; i < m_nc*m_ni*m_nj ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}

void Darray2::writeToFile(char* fileName, int is, int ie, int js, int je)
{
    FILE *extFile = fopen(fileName, "w"); 
    for (int j=js; j<=je; j++)
    {
        for (int i=is; i<=ie; i++)
            fprintf(extFile, " %18.10e", m_data[m_base+m_offi*i+m_offc+m_offj*j]);
        fprintf(extFile,"\n");
    }
    fclose(extFile);
}

