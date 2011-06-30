/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HASH_H_
#define HASH_H_

typedef struct entry entry;
/* Public data for hash element 
	 -- data can either be an [i]nteger or a [p]ointer to something else,
   reverse access to the key from the data can also be useful */
typedef struct {
	union {
    void *p;   
    int i;
	};
	void *key;
  size_t lkey;    /* Length of key */
} data_u;

struct entry { /* Structure of individual entries in list */
  data_u data; /* Data value    */
  unsigned long hashval; /* Cache hash value for use in later expansions */
  entry *next; /* Next item in list */
};

typedef struct {
  unsigned long size, /* Size of active part of hash table */
    frontmask,        /* Mask for front part of hash table */
    fullmask,         /* Mask for both parts of hash table */
    space,            /* Table space allocated at present */
    nelem;            /* Number of elements present */
  void (*freedata)(void *data); /* Function used to free data section */
  entry **tab;        /* Table of entry lists */
  /* Hash function used for this table */
  unsigned long (*hashfunction)(const void *t, const size_t len); 
} hashtab;

/* Create new hash table */
extern hashtab *newhash(void (*freedata)(void *));

/* Remove all of hash table */
extern void freehash(hashtab *table);

/* Look up entry */
extern data_u *lookup(const void *key, size_t lkey, const hashtab *table);

/* Add new entry */
extern data_u *addentry(const void *key, size_t lkey, hashtab *table, int *exists);

/* Determine maximum chain length */
extern int maxchain(const hashtab *table);

/* Make list of data values */
extern unsigned long makelist(const hashtab *table, data_u **list, 
															const unsigned long nlist, int (*maskfun)(data_u *dat));

/* Make list of data values */
extern unsigned long makeplist(const hashtab *table, void **list, 
															const unsigned long nlist, int (*maskfun)(data_u *dat));

#endif /* HASH_H_ */
