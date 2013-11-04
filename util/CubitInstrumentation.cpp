
// File:    CubitInstrumentation.cpp
// Author:  Boyd Tidwell
//
// Purpose: CubitInstrumentation is used to collect and write command and 
//          feature usage by Cubit users. 

#include "CubitInstrumentation.hpp"
#include "CubitMessage.hpp"


//#include <iostream>

CubitInstrumentation* CubitInstrumentation::thisInstance = NULL;

//=============================================================================
CubitInstrumentation* CubitInstrumentation::instance()
{
  if (!thisInstance)
  {
    thisInstance = new CubitInstrumentation;
    if (!thisInstance)
    {
      PRINT_ERROR("Unable to instantiate instrumentation object!\n");
      exit(1);
    }
  }
  return thisInstance;
}

//=============================================================================
void CubitInstrumentation::delete_instance()
{
  if (thisInstance)
    delete thisInstance;
  thisInstance = NULL;
}

//=============================================================================
CubitInstrumentation::CubitInstrumentation()
  : outputState(Unknown), 
    recordTokens(true), validKeywords(false), writeKeywords(true),
    tokenUsageStream(NULL)
{}

//=============================================================================
CubitInstrumentation::~CubitInstrumentation()
{
    // close the token file
  if (tokenUsageStream) 
  {
    close_section();
    tokenUsageStream->close();
    delete tokenUsageStream;
  }
  thisInstance = NULL;
}

//=============================================================================
void CubitInstrumentation::write_context_token(int token)
{
  if (check_token_log())
  {
    if (outputState != GUI)
    {
      outputState = open_section(GUI);
    }
    else
    {
      *tokenUsageStream << ",";
    }
    *tokenUsageStream << token;
    tokenUsageStream->flush();
  }
}

//=============================================================================
void CubitInstrumentation::write_all_words()
{
  if (tokenUsageStream && !keywordMap.empty()) {
    if (outputState != Keyword)
    {
      outputState = open_section(Keyword);
    }

    std::map<std::string, int>::iterator it = keywordMap.begin();
    std::map<std::string, int>::iterator end = keywordMap.end();
    *tokenUsageStream << keywordMap.size() << std::endl;
    *tokenUsageStream << (*it).first;
    ++it;
    for ( ; it != end; it++)
    {
      *tokenUsageStream << "," << (*it).first;
    }
    close_section();
    tokenUsageStream->flush();
  }
}

//=============================================================================
void CubitInstrumentation::load_all_words(const char* const allwords[])
{
  assert(keywordMap.size() == 0);

  int index = 0;
  while (allwords[index] != NULL) {
    keywordMap[allwords[index]] = index + 1;
    ++index;
  }

  write_all_words();
  validKeywords = true;
}

//=============================================================================
int CubitInstrumentation::lookup_keyword(const char* keyword)
{
  if (keyword == NULL)
    return -1;

  std::map<std::string, int>::iterator it = keywordMap.find(keyword);
  if (it != keywordMap.end())
    return (*it).second;

  return -1;
}

//=============================================================================
void CubitInstrumentation::write_keywords(std::vector<CubitString> keywords)
{
  if (check_token_log() && !keywords.empty())
  {
    if (outputState != Command)
    {
      outputState = open_section(Command);
    }
    else
    {
      *tokenUsageStream << ",";
    }
    *tokenUsageStream << lookup_keyword(keywords[0].c_str());
    for (size_t i = 1; i < keywords.size(); i++)
    {
      int debug = lookup_keyword(keywords[i].c_str());
      *tokenUsageStream << ":" << debug;
    }
    tokenUsageStream->flush();
  }
}

//=============================================================================
void CubitInstrumentation::close_section()
{
  if (outputState == Command)
    *tokenUsageStream << std::endl << "</command>" << std::endl;
  else if (outputState == GUI)
    *tokenUsageStream << std::endl << "</gui>" << std::endl;
  else if (outputState == Keyword)
    *tokenUsageStream << std::endl << "</keyword>" << std::endl;
 
  outputState = Unknown;
}

//=============================================================================
CubitInstrumentation::SectionState
CubitInstrumentation::open_section(SectionState state)
{
  close_section();

  if (state == Command)
    *tokenUsageStream << "<command>" << std::endl;
  else if (state == GUI)
    *tokenUsageStream << "<gui>" << std::endl;
  else if (state == Keyword)
    *tokenUsageStream << "<keyword>" << std::endl;
  else
    return Unknown;

  return state;
}

//=============================================================================
bool CubitInstrumentation::check_token_log()
{
  const char* token_file = getenv("_CUBIT_USAGE_TOKEN_FILE");
  if (token_file)
  {
    if (!tokenUsageStream)
    {
      //const char* token_file = getenv("_CUBIT_USAGE_TOKEN_FILE");
      //if(!token_file)
      //  token_file = "usage_tokens.log";
      tokenUsageStream = new std::ofstream(token_file, std::ofstream::out);
      if (!tokenUsageStream)
      {
        std::string file_name = getenv("_CUBIT_USAGE_TOKEN_DIR");
        file_name += std::tmpnam(NULL);
        file_name += "log";
        tokenUsageStream = new std::ofstream(file_name.c_str(), std::ofstream::out);
        if (!tokenUsageStream)
        {
          PRINT_ERROR("Failed to open token usage file!\n");
          PRINT_INFO("\tNo usage token will be recorded.\n");
          return false;
        }
      }
    }

      // write keywords if valid and write flags set
    if (validKeywords && writeKeywords)
    {
      write_all_words();
      writeKeywords = false;
    }
    return true; 
  }
  else
  {
    return false;
  }
}
