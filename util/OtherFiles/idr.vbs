option explicit

private dicFirst

Sub GenerateCubitVersionDotH
'DESCRIPTION: Creates CubitVersion.h with the current date and time
  
  dim datNow, strAmPm, intHour
  dim fso, filDotH, filVersion
  
  ' Make a temp file to write to
  set fso = CreateObject("Scripting.fileSystemObject")
  set filDotH = fso.CreateTextFile("CubitVersion.h", true, false)
  
  ' get the date and time
  datNow = Now
  
  'Fix it up
  intHour = Hour(datNow)
  If intHour > 11 then
    strAmPm = "PM"
    intHour = intHour - 12
  else
    strAmPm = "AM"
  End If
  If intHour = 0 Then intHour = 12
  
  'Write it to the file
  filDotH.WriteLine "#define CUBIT_DATE " & chr(34) & _
     Right("0" & Month(datNow), 2) & "/" & _
     Right("0" & Day(datNow), 2) & "/" & _
     Year(datNow) & " " & _
     Right("0" & intHour, 2) & ":" & _
     Right("0" & Minute(datNow), 2) & " " & strAmPm & chr(34)
  
  'We will get the non-date related version number from CubitVersion.tmp
  Set filVersion = fso.OpenTextFile("nt_gui_source\CubitVersion.tmp", 1)
  
  'Write out the cubit version
  filDotH.writeline "#define CUBIT_VERSION " & chr(34) & filVersion.ReadLine & "-" & _
     Right(Year(datNow), 2) & _
     Right("00" & DatePart("y", datNow), 3) & chr(34)
  
  filVersion.close
  filDotH.close
  
End Sub

private Function GetCommandTableNames()
'DESCRIPTION: Gets the files listed in COMMAND_TABLES in the makefile.
'             Returns the names as a string array

  dim fso, filMakefile
  dim str, intCount, boolContinue
  dim strFileArray(25)
  
  intCount = 0
  
  ' Get the list of command tables from the makefile
  set fso = CreateObject("Scripting.FileSystemObject")
  set filMakefile = fso.OpenTextFile("makefile", 1)
  str = ""
  while left(str, 14) <> "COMMAND_TABLES"
    str = filMakefile.readline
  wend
  
  ' str now holds COMMAND_TABLES = \
  ' make sure it doesn't have a file name too
  boolContinue = (right(str, 1) = "\")
  str = trim(mid(str, instr(str, "=")+1))
  If str <> "\" and str <> "" then
    strFileArray(0) = trim(left(str, len(str)-1))
    intCount = 1
  End If
  
  while boolContinue = true
    'Get the next line
    str = filmakefile.readline
    
    'Get rid of the trailing "\"
    If right(str,1) = "\" then
      str = left(str, len(str)-1)
      boolContinue = true
    else
      boolContinue = false
    End If
    
    ' Get rid of any leading tabs or spaces
    while left(str,1) = " " or left(str, 1) = "	"
      str = mid(str,2)
    wend
    
    ' Get rid of any trailing tabs or spaces
    while right(str,1) = " " or right(str, 1) = "	"
      str = left(str, len(str) - 1)
    wend

    ' If this is a variable name instead of a file name,
    ' ignore it.  This isn't the safest thing to do, but
    ' it works for now.
    If left(str,1) <> "$" then
      ' Save the name of the file
      strFileArray(intCount) = str
      intCount = intCount + 1
    End If
    
  wend

  GetCommandTableNames = strFileArray
  
End Function

private sub CreateDotH(strFileName, boolKeywordsOnly)
  dim fso, strFileNames', txfTmpFile
  dim intFileIndex, filcurFile
  dim str, intPos1, intPos2
  

  ' Get the files to look in
  strFileNames = GetCommandTableNames
  intFileIndex = 0
  While intFileIndex < UBound(strFileNames) And Not IsEmpty(strFileNames(intFileIndex))
    intFileIndex = intFileIndex + 1
  Wend
  WScript.Echo " " & intFileIndex & " Files to Process"
  
  ' Open a temp file
  set fso = CreateObject("Scripting.FileSystemObject")
  'set txfTmpFile = fso.CreateTextFile(strFileName, true, false)

  'Initialize the sorting tables
  InitializeWordArrays
  
  'Go through each line of each file
  intFileIndex = 0
  While Not IsEmpty(strFileNames(intFileIndex))
    WScript.Echo "  Processing File " & intFileIndex + 1 & ": " & strFileNames(intFileIndex)
    set filCurFile = fso.OpenTextFile(strFileNames(intFileIndex), 1, 0)
    while filCurFile.atendofstream = false
      str = filCurFile.readline
      If boolKeywordsOnly then
        do while len(str) > 0
          ' Look for an opening bracket
          intPos1 = instr(str, "{")
          If intPos1 = 0 then exit do
          str = mid(str, intPos1 + 1)

          ' Look for a beginning quote
          intPos1 = instr(str, chr(34))
          If intPos1 = 0 then exit do
          
          'Find the end quote
          intPos2 = instr(intPos1+1, str, chr(34))
          ' Make sure it isn't an empty string or a string that doesn't start with
          ' a letter
          If (intPos2 > intPos1 + 1) and _
             (strcomp(mid(str, intPos1+1, 1), "A", 1) >= 0) and _
             (strcomp(mid(str, intPos1+1, 1), "Z", 1) <= 0) then
            ' Write it to a file
            'txfTmpFile.WriteLine (mid(str, intPos1, intPos2 - intPos1 + 1) & ",")
            AddWordToArrays mid(str, intPos1+1, intPos2-intPos1-1)
          End If
          ' Set the string to the portion we haven't read yet
          str = mid(str, intPos2 + 1)
        loop
      else
        do while len(str) > 0
          ' Look for a beginning quote
          intPos1 = instr(str, chr(34))
          If intPos1 = 0 then exit do

          ' See If there was a #include before the quotes
          intPos2 = instr(str, "#include")
          If intPos2 <> 0 and intPos2 < intPos1 then
            'Skip the first set of quotes after the #include
            intPos2 = instr(intPos1+1, str, chr(34))
            str = mid(str, intPos2 + 1)
          else
            ' There was no #include, so use the string
            intPos2 = instr(intPos1+1, str, chr(34))
            ' Make sure it isn't an empty string or a string that doesn't start with
            ' a letter
            If (intPos2 > intPos1 + 1) and _
               (strcomp(mid(str, intPos1+1, 1), "A", 1) >= 0) and _
               (strcomp(mid(str, intPos1+1, 1), "Z", 1) <= 0) then
              ' Write it to a file
              'txfTmpFile.WriteLine (mid(str, intPos1, intPos2 - intPos1 + 1) & ",")
              AddWordToArrays mid(str, intPos1+1, intPos2-intPos1-1)
            End If
            ' Set the string to the portion we haven't read yet
            str = mid(str, intPos2 + 1)
          End If
        loop
      End If
    Wend
    intFileIndex = intFileIndex + 1
    filCurFile.close
  Wend
  
  'Now write the file
  WriteFileFromArrays strFileName

End Sub

private Sub InitializeWordArrays
  set dicFirst = CreateObject("Scripting.Dictionary")
  dicFirst.add NULL, ""
End Sub

private sub AddWordToArrays(strWord)
  dim dicNew, dicCur, dicPrev
  dim vals, keyArray, intCompVal

  'Printtooutputwindow("Adding word " & strWord)

  'Make a new array entry for this word
  set dicNew = CreateObject("Scripting.Dictionary")
  dicNew.add NULL, strWord

  'Set dicCur to the first word, and dicPrev to the head
  set dicPrev = dicFirst
  keyArray = dicPrev.keys
  If IsNull(keyArray(0)) = false then
    set dicCur = keyArray(0)
  else
    dicCur = NULL
  End If

  'Find where the new word fits in the list
  do while IsNull(dicCur) <> true
    vals = dicCur.items
    intCompVal = StrComp(strWord, vals(0), 0)
    If intCompVal = 0 then
      'Skip this word...it is a repeat
      set dicNew = Nothing
      exit sub
    elseif intCompVal = -1 then
      'dicNew comes before dicCur
      exit do
    End If

    'If we're still here, we need to keep looking
    keyArray = dicCur.keys
    set dicPrev = dicCur
    If IsNull(keyArray(0)) = false then
      set dicCur = keyArray(0)
    else
      dicCur = NULL
    End If
  loop
  dicNew.key(NULL) = dicCur
  dicPrev.key(dicCur) = dicNew
  
end sub

private sub WriteFileFromArrays(strFileName)
  dim txfOut, vals, keyArray, dicCur, fso

  'Open the file for writing
  set fso = CreateObject("Scripting.FileSystemObject")
  set txfOut = fso.CreateTextFile(strFileName, true, false)

  'Start at the beginning of the arrays
  keyArray = dicFirst.keys
  If IsNull(keyArray(0)) = true then
    dicCur = NULL
  else
    set dicCur = keyArray(0)
  End If

  do while IsNull(dicCur) <> true
    vals = dicCur.items
    keyArray = dicCur.keys
    txfOut.WriteLine chr(34) & vals(0) & chr(34) & ", "
    If IsNull(keyArray(0)) = true then
      dicCur = NULL
    else
      set dicCur = keyArray(0)
    End If
  loop

  txfOut.close

End Sub

Sub TouchCubitVersion()
'DESCRIPTION: Resets the modification date of CubitVersion.tmp

  Dim fso, txf1, txf2
  
  Set fso = CreateObject("Scripting.FileSystemObject")
  
  'Make sure we're compiling CUBIT
  If ActiveProject.Name = "Cubit" then

    'Change directories to nt_gui_source
    Application.CurrentDirectory = _
      fso.GetParentFolderName(_
        fso.GetParentFolderName(ActiveProject.FullName))

    'Make sure CubitVersion.tmp exists
    If  fso.FileExists("CubitVersion.tmp") Then

      'Create new file, copy old to new
      Set txf1 = fso.OpenTextFile("CubitVersion.tmp", 1)
      Set txf2 = fso.CreateTextFile("CubitVersionTmp.tmp", True)
      txf2.Write txf1.ReadAll
      txf1.Close
      txf2.Close

      'Delete the old file, Rename the new file
      fso.DeleteFile "CubitVersion.tmp"
      fso.MoveFile "CubitVersionTmp.tmp", "CubitVersion.tmp"
      
    End If
  End If
End Sub


Sub GenerateAllwordsDotH
  'DESCRIPTION: Creates allwords.h from the command tables.
  WScript.Echo "Generating allwords.h..."
  CreateDotH "allwords.h", false
  WScript.Echo "Generation of allwords.h is complete"
End Sub

Sub GenerateKeywordsDotH
  'DESCRIPTION: Creates keywords.h from the command tables.
  WScript.Echo "Generating keywords.h..."
  CreateDotH "keywords.h", true
  WScript.Echo "Generation of keywords.h is complete"
End Sub

dim objArgs

set objArgs = WScript.Arguments

If objArgs.Count < 1 Then
  WScript.Echo "  Error: At least one command-line argument required"
End If

Select Case objArgs(0)
  case "/GenerateAllwordsDotH"
    GenerateAllwordsDotH
  case "/GenerateKeywordsDotH"
    GenerateKeywordsDotH
End Select
