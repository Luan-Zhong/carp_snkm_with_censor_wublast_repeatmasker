import sys
import re

BlastFileName = sys.argv[1]

def StartHit():
	global QueryFromPos
	global QueryToPos
	global SubjectFromPos
	global SubjectToPos
	global SubjectId
	global PctId
	global Evalue

	EndHit()

	QueryFromPos = 0
	QueryToPos = 0
	SubjectFromPos = 0
	SubjectToPos = 0
	Evalue = 0.0
	PctId = -999.0

def IdFromAnnot(Annot):
	Fields = Annot.split("|")
	if len(Fields) < 2:
		Fields = Annot.split(" ")
		return Fields[0]
	return Fields[0] + "|" + Fields[1]

def FreeTextFromAnnot(Annot):
	Pos = Annot.find(" ")
	if Pos < 0:
		return ""
	Text = Annot[Pos+1:]
	if len(Text) > 32:
		Text = Text[:32] + "..."
	return Text

def SubjectLengthLine(Line):
	global SubjectLength
	Fields = Line.split("=")
	if len(Fields) != 2:
		print >> sys.stderr, "** Expected Length=, got: ", Line
		sys.exit(1)
	SubjectLength = int(Fields[1].strip().replace(",", ""))

def QueryLengthLine(Line):
	global QueryLength
	Fields = Line.split("(")
	if len(Fields) != 2:
		print >> sys.stderr, "** Expected Length=, got: ", Line
		sys.exit(1)
	Fields2 = Fields[1].split(" ")
	strQueryLength = Fields2[0]
#	print Line, Fields, "'" + strQueryLength + "'"
	QueryLength = int(strQueryLength.replace(",", ""))

def StrandLine(Line):
	global Rev
	Line = Line.strip()
	if Line == "Strand = Plus / Plus":
		Rev = 0
	elif Line == "Strand = Plus / Minus":
		Rev = 1
	else:
		print >> sys.stderr, "** Unexpected strand:", Line
		sys.exit(1)

def SubjectIdLine(Line):
	global SubjectId
	global SubjectLength
	global SubjectFreeText
	EndHit()
	Annot = Line[1:].strip()
	SubjectId = IdFromAnnot(Annot)
	SubjectFreeText = FreeTextFromAnnot(Annot)
	SubjectLength = -999

def QueryIdLine(Line):
	global QueryId
	EndHit()
	QueryId = IdFromAnnot(Line[7:].strip())

def ScoreLine(Line):
	global Evalue
	StartHit()
	Fields = Line.split("=")
	if len(Fields) != 4:
		print "** ERROR ** ScoreLine ", Line, " len(Fields)=", len(Fields)
		sys.exit(1)
	Fields2 = Fields[2].split(",")
	e = Fields2[0].strip()
	if e[0] == "e":
		e = "1" + e
	Evalue = float(e)

def QueryLine(Line):
	global QueryFromPos
	global QueryToPos

	s = Line[7:-1]
	Fields = s.split()
	n = len(Fields)
	if n < 3:
		print "** ERROR ** < 3 fields s='" + s + "'"
		sys.exit(1)

	FromPos = int(Fields[0])
	ToPos = int(Fields[n-1])

	if QueryFromPos == 0:
		QueryFromPos = FromPos
	QueryToPos = ToPos

def SubjectLine(Line):
	global SubjectFromPos
	global SubjectToPos

	s = Line[7:-1]
	Fields = s.split()
	n = len(Fields)
	if n < 3:
		print "** ERROR ** < 3 fields s='" + s + "'"
		sys.exit(1)

	FromPos = int(Fields[0])
	ToPos = int(Fields[n-1])

	if SubjectFromPos == 0:
		SubjectFromPos = FromPos
	SubjectToPos = ToPos
	
def EndHit():
	global QueryFromPos
	global QueryToPos
	global SubjectFromPos
	global SubjectToPos
	global QueryId
	global SubjectId
	global SubjectFreeText
	global QueryLength
	global SubjectLength
	global Rev
	global Evalue

	if QueryFromPos == 0:
		return

	if Rev:
		Strand = "-"
		Tmp = SubjectFromPos
		SubjectFromPos = SubjectToPos
		SubjectToPos = Tmp
	else:
		Strand = "+"
		
	if QueryFromPos > QueryToPos:
		if Rev == 1:
			Strand = "+"
		else:
			Strand = "-"
		Tmp = QueryFromPos
		QueryFromPos = QueryToPos
		QueryToPos = Tmp

	if SubjectFromPos > SubjectToPos:
		Tmp = SubjectFromPos
		SubjectFromPos = SubjectToPos
		SubjectToPos = Tmp
	
	print "%s\tblast\thit\t%d\t%d\t%.3g\t%s\t.\tTarget %s %d %d; QueryLength %d; TargetLength %d; Annot \"%s\"" % \
	  ( QueryId, QueryFromPos, QueryToPos, Evalue, Strand, SubjectId, SubjectFromPos, SubjectToPos, QueryLength, SubjectLength, SubjectFreeText)

	QueryFromPos = 0

def IdentitiesLine(Line):
	global PctId

	Fields1 = Line[14:].split(" ")
	Fields2 = Fields1[0].split("/")
	if len(Fields2) != 2:
		print "** ERROR ** Idents line (Fields2) ", Line
		sys.exit(1)
	x = int(Fields2[0])
	y = int(Fields2[1])
	PctId = (x*100.0)/y

# End PSI-BLAST iterations
def RoundLine(Line):
	return

QueryId = ""
SubjectId = ""
SubjectFreeText = ""
QueryFromPos = 0
QueryLength = 0

StartHit()

Rev = 0
File = open(BlastFileName)
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break

	if Line.startswith(">"):
		SubjectIdLine(Line)
	elif Line.startswith("Query="):
		QueryIdLine(Line)
	elif Line.startswith("Query:"):
		QueryLine(Line)
	elif Line.startswith("Sbjct:"):
		SubjectLine(Line)
	elif Line.startswith(" Identities ="):
		IdentitiesLine(Line)
	elif Line.startswith(" Score = "):
		ScoreLine(Line)
	elif Line.startswith("Results from round"):
		RoundLine(Line)
#	elif Line.startswith("        Length ="):
	elif re.search("Length =", Line) != None:
		SubjectLengthLine(Line)
#	elif Line.strip().endswith("letters)"):
	elif re.search(" letters;", Line) != None:
		QueryLengthLine(Line);
	elif Line.startswith(" Strand ="):
		StrandLine(Line)

EndHit()
