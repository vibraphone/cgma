#include "stdafx.h"
#include <assert.h>

HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, double *dArray)
{
    HRESULT hr = NOERROR;
    double *dData;

    // make sure we can copy the full array
    long lUBound;
    long lLBound;
    SAFEARRAY* psa = V_ARRAY(&v);
    hr = SafeArrayGetUBound(psa, 1, &lUBound);
    hr = SafeArrayGetLBound(psa, 1, &lLBound);
    if ((lUBound - lLBound + 1) > lArraySize )
    {
        assert(false);
        return E_FAIL;
    }

    // access the safearray data
    hr = SafeArrayAccessData(psa, (void**)&dData);
    if (FAILED(hr)) return hr;

    // copy the data from the safearray to the input array
	for (int i=lLBound, j=0; i<lUBound+1; i++, j++) dArray[j] = dData[i];

    // Unaccess the component SafeArray
    hr = SafeArrayUnaccessData(psa);  

    return S_OK;
}

HRESULT GetVARIANTArraySize(VARIANT &v, long &lArraySize)
{
    HRESULT hr = NOERROR;
    SAFEARRAY* psa = V_ARRAY(&v);

    long lUBound;
    long lLBound;
    hr = SafeArrayGetUBound(psa, 1, &lUBound);
    hr = SafeArrayGetLBound(psa, 1, &lLBound);

    lArraySize = lUBound - lLBound + 1;
    return hr;
}

HRESULT ExtractVARIANTArrayData(VARIANT &v, long lArraySize, LPDISPATCH *pDispArray)
{
    HRESULT hr = NOERROR;
    LPDISPATCH *pDispData;

    // make sure we can copy the full array
    long lUBound;
    long lLBound;
    SAFEARRAY* psa = V_ARRAY(&v);
    hr = SafeArrayGetUBound(psa, 1, &lUBound);
    hr = SafeArrayGetLBound(psa, 1, &lLBound);
    if ((lUBound - lLBound + 1) > lArraySize )
    {
        assert(false);
        return E_FAIL;
    }

    // access the safearray data
    hr = SafeArrayAccessData(psa, (void**)&pDispData);
    if (FAILED(hr)) return hr;

    // copy the data from the safearray to the input array
	for (int i=lLBound, j=0; i<lUBound+1; i++, j++) pDispArray[j] = pDispData[i];

    // Unaccess the component SafeArray
    hr = SafeArrayUnaccessData(psa);  

    return S_OK;
}