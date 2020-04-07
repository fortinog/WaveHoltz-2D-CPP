#include <iostream>

#include "Darray.h"

using namespace std;

// Default value of array ordering
bool Darray::m_corder = false;

// Constructors
//-----------------------------------------------------------------------
Darray::Darray( int nc, int ibeg, int iend)
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
Darray::Darray(int ibeg, int iend)
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
Darray::Darray()
{
    m_nc = m_ib = m_ie = 0;
    m_data = NULL;
}

//
void Darray::define( int nc, int ibeg, int iend)
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

void Darray::define(int ibeg, int iend)
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
void Darray::set_value(double scalar)
{
    for( size_t i = 0 ; i < m_npts ; i++ )
        m_data[i] = scalar;
}

//-----------------------------------------------------------------------
void Darray::define_offsets()
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
void Darray::copy( const Darray& u )
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


