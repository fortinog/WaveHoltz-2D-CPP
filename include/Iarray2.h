// A simple two dimensional array class that can use arbitrary index range
#ifndef IARRAY2
#define IARRAY2

class Iarray2
{
  public:
    Iarray2( int nc, int x_s, int x_e, int y_s, int y_e);
    Iarray2( int x_s, int x_e, int y_s, int y_e);
    Iarray2();
    ~Iarray2() {if( m_data != 0 ) delete[] m_data;}
    void define( int nc, int x_s, int x_e, int y_s, int y_e);
    void define( int x_s, int x_e, int y_s, int y_e);
    inline bool in_range( int c, int i, int j)
    {return 1 <= c && c <= m_nc && m_s <= i && i <= m_e && n_s <= j && j <= n_e;}

// overload parenthesis operator 
    inline int& operator()( int c, int i, int j){
            // Turn on array bounds check 
#ifdef BZ_DEBUG  
        try
        {  
            if (!in_range(c,i,j)) throw 10;
        }
        catch(int e) 
        {
            std::cout <<
                "Error Index (c,i,j) = (" << c << "," <<
                i << ") not in range 1<= c <= " <<
                m_nc << " " << m_s << " <= i <= " << m_e
                << " " << n_s << " <= j <= " << n_e;
        }
#endif
       return m_data[m_base+off_c*c+off_i*i +off_j*j];
    }
    inline int& operator()( int i, int j)
    {
#ifdef BZ_DEBUG
        try{
            if (!in_range(1,i,j)) throw 10;
        }
        catch(int e)
        {
            std::cout << "Error Index (c,i,j,k) = (" << 1
                      << "," << i << ") not in range 1<= c <= "
                      << m_nc << " "
                      << m_s << " <= i <= " << m_e;
        }
#endif
    return m_data[m_base+off_c+off_i*i +off_j*j];
    }
    inline bool is_defined(){return m_data != NULL;}
    static bool m_corder;

    // Beginning and ending index along each dimension
    // Note *_s means start, *_e means end
    int m_s, m_e, n_s, n_e;
    ssize_t m_base;
    size_t off_i, off_j, off_c, npts;
    void define_offsets();
    void set_value( int scalar );
    void copy( const Iarray2& u );
    int m_nc, M, N;
  private:
   int* m_data;
};
#endif