#include "SWRefKey.h"
#include "CubitDefines.h"
#include <memory.h>

SWRefKey::SWRefKey()
{
    m_nBytes = 0;
    m_pKeyData = NULL;
}

SWRefKey::SWRefKey(const SWRefKey &other)
{
    m_nBytes = other.m_nBytes;;
    m_pKeyData = new unsigned char[m_nBytes];
    assert(m_pKeyData);
    memcpy(m_pKeyData, other.m_pKeyData, m_nBytes);
}

SWRefKey::~SWRefKey()
{
    if (m_pKeyData)
        delete [] m_pKeyData;
    m_pKeyData = NULL;
    m_nBytes = 0;
}

SWRefKey &SWRefKey::operator= (const SWRefKey &other)
{
    // clean up old data
    if (m_pKeyData)
    {
        delete [] m_pKeyData;
        m_pKeyData = NULL;
        m_nBytes = 0;
    }

    // copy data
    m_nBytes = other.m_nBytes;;
    m_pKeyData = new unsigned char[m_nBytes];
    assert(m_pKeyData);
    memcpy(m_pKeyData, other.m_pKeyData, m_nBytes);
    return *this;
}

unsigned long SWRefKey::putData(unsigned long nBytes, unsigned char* pData)
{
    assert(nBytes);
    if (!nBytes) return 0;
    m_nBytes = nBytes;
    m_pKeyData = new unsigned char[nBytes];
    if (m_pKeyData)
    {
        memcpy(m_pKeyData, pData, nBytes);
        return m_nBytes;
    }
    else
    {
        return 0;
    }
}

unsigned long SWRefKey::size()
{
    return m_nBytes;
}

unsigned char *SWRefKey::keyData()
{
    return m_pKeyData;
}

