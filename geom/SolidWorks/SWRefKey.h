#ifndef __SWREFKEY_H__
#define __SWREFKEY_H__

class SWRefKey
{
public:
    SWRefKey();
    SWRefKey(const SWRefKey &other);
    ~SWRefKey();

    SWRefKey &operator= (const SWRefKey &other);

    unsigned long putData(unsigned long nBytes, unsigned char * pData);

    unsigned long size();
    unsigned char *keyData();

private:
	unsigned long m_nBytes;
	unsigned char * m_pKeyData;
};
#endif