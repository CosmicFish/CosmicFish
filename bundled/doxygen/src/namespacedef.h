/******************************************************************************
 *
 * 
 *
 * Copyright (C) 1997-2014 by Dimitri van Heesch.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation under the terms of the GNU General Public License is hereby 
 * granted. No representations are made about the suitability of this software 
 * for any purpose. It is provided "as is" without express or implied warranty.
 * See the GNU General Public License for more details.
 *
 * Documents produced by Doxygen are derivative works derived from the
 * input used in their production; they are not affected by this license.
 *
 */

#ifndef NAMESPACEDEF_H
#define NAMESPACEDEF_H

#include <qstrlist.h>
#include <qdict.h>
#include "sortdict.h"
#include "definition.h"
#include "filedef.h"

class MemberList;
class ClassDef;
class ClassList;
class OutputList;
class ClassSDict;
class MemberDef;
class NamespaceList;
class MemberGroupSDict;
class NamespaceSDict;
class FTextStream;

/** A model of a namespace symbol. */
class NamespaceDef : public Definition
{
  public:
    NamespaceDef(const char *defFileName,int defLine,int defColumn,
                 const char *name,const char *ref=0,
                 const char *refFile=0,const char*type=0,
                 bool isPublished=false);
   ~NamespaceDef();
    DefType definitionType() const { return TypeNamespace; }
    QCString getOutputFileBase() const;
    QCString anchor() const { return QCString(); }
    void insertUsedFile(FileDef *fd);
    
    void writeDocumentation(OutputList &ol);
    void writeMemberPages(OutputList &ol);
    void writeQuickMemberLinks(OutputList &ol,MemberDef *currentMd) const;
    void writeTagFile(FTextStream &);

    void insertClass(ClassDef *cd);
    void insertNamespace(NamespaceDef *nd);
    void insertMember(MemberDef *md);

    void computeAnchors();
    int countMembers();
    void addUsingDirective(NamespaceDef *nd);
    NamespaceSDict *getUsedNamespaces() const;
    void addUsingDeclaration(Definition *def);
    SDict<Definition> *getUsedClasses() const { return usingDeclList; }
    void combineUsingRelations();
    QCString displayName(bool=TRUE) const;
    QCString localName() const;

    bool isConstantGroup() const { return CONSTANT_GROUP == m_type; }
    bool isModule()        const { return MODULE == m_type; }
    bool isLibrary() const { return LIBRARY == m_type; }

    bool isLinkableInProject() const;
    bool isLinkable() const;
    bool hasDetailedDescription() const;
    void addMembersToMemberGroup();
    void distributeMemberGroupDocumentation();
    void findSectionsInDocumentation();
    void sortMemberLists();

    virtual Definition *findInnerCompound(const char *name);
    void addInnerCompound(Definition *d);
    void addListReferences();

    bool subGrouping() const { return m_subGrouping; }
    
    MemberList *getMemberList(MemberListType lt) const;
    const QList<MemberList> &getMemberLists() const { return m_memberLists; }
    MemberDef    *getMemberByName(const QCString &) const;

    /*! Returns the user defined member groups */
    MemberGroupSDict *getMemberGroupSDict() const { return memberGroupSDict; }

    /*! Returns the classes contained in this namespace */
    ClassSDict *getClassSDict() const { return classSDict; }

    /*! Returns the namespaces contained in this namespace */
    NamespaceSDict *getNamespaceSDict() const { return namespaceSDict; }

    QCString title() const;
    QCString compoundTypeString() const;

    bool visited;

  private:
    MemberList *createMemberList(MemberListType lt);
    void addMemberToList(MemberListType lt,MemberDef *md);
    void writeMemberDeclarations(OutputList &ol,MemberListType lt,const QCString &title);
    void writeMemberDocumentation(OutputList &ol,MemberListType lt,const QCString &title);
    void writeDetailedDescription(OutputList &ol,const QCString &title);
    void writeBriefDescription(OutputList &ol);
    void startMemberDeclarations(OutputList &ol);
    void endMemberDeclarations(OutputList &ol);
    void writeClassDeclarations(OutputList &ol,const QCString &title);
    void writeInlineClasses(OutputList &ol);
    void writeNamespaceDeclarations(OutputList &ol,const QCString &title,
            bool isConstantGroup=false);
    void writeMemberGroups(OutputList &ol);
    void writeAuthorSection(OutputList &ol);
    void startMemberDocumentation(OutputList &ol);
    void endMemberDocumentation(OutputList &ol);
    void writeSummaryLinks(OutputList &ol);
    void addNamespaceAttributes(OutputList &ol);

    QCString              fileName;
    FileList              files;

    NamespaceSDict       *usingDirList;
    SDict<Definition>    *usingDeclList;
    SDict<Definition>    *m_innerCompounds;

    MemberSDict          *m_allMembersDict;
    QList<MemberList>     m_memberLists;
    MemberGroupSDict     *memberGroupSDict;
    ClassSDict           *classSDict;
    NamespaceSDict       *namespaceSDict;
    bool                  m_subGrouping;
    enum { NAMESPACE, MODULE, CONSTANT_GROUP, LIBRARY } m_type;
    bool m_isPublished;
};

/** A list of NamespaceDef objects. */
class NamespaceList : public QList<NamespaceDef>
{
  public:
   ~NamespaceList() {}
    int compareValues(const NamespaceDef *nd1,const NamespaceDef *nd2) const
    {
      return qstricmp(nd1->name(), nd2->name());
    }
};

/** An iterator for NamespaceDef objects in a NamespaceList. */
class NamespaceListIterator : public QListIterator<NamespaceDef>
{
  public:
    NamespaceListIterator(const NamespaceList &l) : 
      QListIterator<NamespaceDef>(l) {}
};

/** An unsorted dictionary of NamespaceDef objects. */
class NamespaceDict : public QDict<NamespaceDef>
{
  public:
    NamespaceDict(int size) : QDict<NamespaceDef>(size) {}
   ~NamespaceDict() {}
};

/** A sorted dictionary of NamespaceDef objects. */
class NamespaceSDict : public SDict<NamespaceDef>
{
  public:
    NamespaceSDict(int size=17) : SDict<NamespaceDef>(size) {}
   ~NamespaceSDict() {}
    void writeDeclaration(OutputList &ol,const char *title,
            bool isConstantGroup=false, bool localName=FALSE);
    bool declVisible() const;
  private:
    int compareValues(const NamespaceDef *item1,const NamespaceDef *item2) const
    {
      return qstricmp(item1->name(),item2->name());
    }
};



#endif
