#include <iostream>
#include <sys/types.h>
#include "Iarray2.h"

using namespace std;

// Default value of array ordering
bool Iarray2::m_corder = true;

// Constructors
//-----------------------------------------------------------------------
Iarray2::Iarray2( int nc, int x_s, int x_e, int y_s, int y_e)
{
    m_nc = nc;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;
    M    = m_e-m_s+1;
    N    = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new int[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}
// Single variable array
Iarray2::Iarray2(int x_s, int x_e, int y_s, int y_e)
{
    m_nc = 1;
    m_s  = x_s;
    m_e  = x_e;
    n_s  = y_s;
    n_e  = y_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;
    if( m_nc*M*N > 0 )
        m_data = new int[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}

//-----------------------------------------------------------------------

//Default constructor
Iarray2::Iarray2()
{
    m_nc = m_s = m_e = n_s = n_e = 0;
    m_data = NULL;
}

//
void Iarray2::define( int nc, int x_s, int x_e, int y_s, int y_e)
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
        m_data = new int[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}

void Iarray2::define(int x_s, int x_e, int y_s, int y_e)
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
        m_data = new int[m_nc*M*N];
    else
        m_data = NULL;
    define_offsets();
}


//-----------------------------------------------------------------------
void Iarray2::set_value(int scalar)
{
    for( size_t i = 0 ; i < npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Iarray2::define_offsets()
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
void Iarray2::copy( const Iarray2& u )
{
    if( m_data != NULL )
        delete[] m_data;
    
    m_nc = u.m_nc;
    m_s = u.m_s;
    m_e = u.m_e;
    n_s = u.n_s;
    n_e = u.n_e;

    M = m_e-m_s+1;
    N = n_e-n_s+1;

    if( m_nc*M*N > 0 )
    {
        m_data = new int[m_nc*M*N];
        for( int i = 0 ; i < m_nc*M*N ; i++ )
            m_data[i] = u.m_data[i];
    }
    else
        m_data = NULL;
    define_offsets();
}
