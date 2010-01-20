//-------------------------------------------------------------------------
// Filename      : OrderedMap.hpp
//
// Purpose       : An std::map that keeps order that things were inserted
//                 into it.
//
// Creator       : Matt Staten
//
// Creation Date : 01/27/2006
//
// Owner         : Matt Staten
//-------------------------------------------------------------------------

#ifndef OrderedMap_HPP
#define OrderedMap_HPP

#include <map>
#include "DLIList.hpp"

template<class X, class Y> class OrderedMap
{
private:
    DLIList<X> mList1;
    std::map<X,Y> mMap;

public:

    OrderedMap( void ) {};

    bool find( X key, Y *data = NULL ) const
    {
        typename std::map<X,Y>::const_iterator iter = mMap.find( key );
        if ( iter == mMap.end() )
        {
            return CUBIT_FALSE;
        }
        if ( data )
        {
            *data = iter->second;
        }
        return CUBIT_TRUE;
    };

    // Return true if the new pair was inserted.
    // Return false if this key is already being used.
    bool insert( X new_key, Y new_data )
    {
        if ( !find( new_key ) )
        {
            mList1.append( new_key );
            mMap.insert( std::map<X,Y>::value_type( new_key, new_data ) );
            return CUBIT_TRUE;
        }
        return CUBIT_FALSE;
    };

    // This can be inefficient
    void erase( X key_to_erase )
    {
        typename std::map<X,Y>::iterator iter = mMap.find( key_to_erase );
        if ( iter != mMap.end() )
        {
            mMap.erase( iter );
            mList1.remove_all_with_value( key_to_erase );
        }
    };

    void erase_iter( int iter_to_erase )
    {
        X key_to_erase = mList1[iter_to_erase];
        typename std::map<X,Y>::iterator iter = mMap.find( key_to_erase );
        if ( iter != mMap.end() )
        {
            mMap.erase( iter );
            mList1.remove_all_with_value( key_to_erase );
        }
    };

    X first( int iter ) const
    {
        return mList1[iter];
    }
    Y second( int iter ) const
    {
        Y data;
        bool stat = find( mList1[iter], &data );
        assert(stat);
        return data;
    }

    bool empty( void ) const { return ( mList1.size() == 0 ? CUBIT_TRUE : CUBIT_FALSE ); };
    int size( void ) const { return mList1.size(); };

    void clean_out( void )
    {
        mList1.clean_out();
        mMap.clear();
    };

    OrderedMap<X,Y>& operator=(const OrderedMap<X,Y>&from)
    {
    	clean_out();
    	for ( int i = 0; i < from.size(); i++ )
    	{
    		X key = from.first( i );
    		Y data = from.second( i );
    		this->insert( key, data );
    	}
    	return *this;
    };

    Y operator[]( X key ) const
    {
        Y data;
        bool stat = find( key, &data );
        assert(stat);
        return data;
    };
};

#endif // OrderedMap_HPP
