#include <iostream>

#include "Darray2.h"

using namespace std;

// Default value of array ordering
bool Darray2::m_corder = true;

// Constructors
//-----------------------------------------------------------------------
Darray2::Darray2( int nc, int x_s, int x_e, int y_s, int y_e)
{
    m_nc = nc;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;
    M    = m_e-m_s+1;
    N    = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new double[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
Darray2::Darray2(int x_s, int x_e, int y_s, int y_e)
{
    m_nc = 1;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new double[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
Darray2::Darray2()
{
    m_nc = m_s = m_e = n_s = n_e = 0;
    m_data = NULL;
}

//
void Darray2::define( int nc, int x_s, int x_e, int y_s, int y_e)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = nc;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new double[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}

void Darray2::define(int x_s, int x_e, int y_s, int y_e)
{
    if( m_data != NULL )
        delete[] m_data;

    m_nc = 1;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new double[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}


//-----------------------------------------------------------------------
void Darray2::set_value(double scalar)
{
    for( size_t i = 0 ; i < npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Darray2::define_offsets()
{
    npts = static_cast<size_t>(M)*m_nc*(N);
    if( m_corder )
    {
        off_c = M*N;
        off_i = N;
        off_j = 1;
        m_base = -m_s*off_i - n_s - off_c;
    }
    else
    {
            // (c,i)=c-1+nc*(i-ib)
        // FORTRAN ordering
        off_c = M*N;
        off_i = 1;
        off_j = M;
        m_base = -n_s*off_j - m_s - off_c;
    }
}

//-----------------------------------------------------------------------
void Darray2::copy( const Darray2& u )
{
    // if( m_data != NULL )
    //     delete[] m_data;
    
    m_nc = u.m_nc;
    m_s = u.m_s;
    m_e = u.m_e;
    n_s = u.n_s;
    n_e = u.n_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;

    if( m_nc*M*N > 0 )
    {
        // m_data = new double[m_nc*M*N];
        for( int i = 0 ; i < m_nc*M*N ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}

void Darray2::writeToFile(char* fileName, int is, int ie, int js, int je)
{
    FILE *extFile = fopen(fileName, "w");
    for (int i=is; i<=ie; i++) 
    {
        for (int j=js; j<=je; j++)
            fprintf(extFile, " %18.10e", m_data[m_base+off_i*i+off_c+off_j*j]);
        fprintf(extFile,"\n");
    }
    fclose(extFile);
}