/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hash.h"

/* Implement sorted linear hash table -- automatically expands to
 * get uniform bin filling, requires good hash function as 2^n masking is used.
 *
 * References
 *
 * W. Litwin, Linear hashing: A new tool for file and table addressing, 
 * Proc. 6th Conference on Very Large Databases, pages 212-223, 1980.
 *
 * http://www.concentric.net/~Ttwang/tech/sorthash.htm
 *
 */

/* Function to calculate 4 byte hash index -- must be good for 2^n table */
STATIC unsigned long hashfunction(const void *t, const size_t len);

/* Internal utility functions for hash table */
STATIC entry *lookup_key(const void *key, size_t *lkey, const hashtab *table, 
												 unsigned long *hk);
STATIC void extend(hashtab *table);
STATIC unsigned long getbin(unsigned long hk, const hashtab *table);


hashtab *newhash(void (*freedata)(void *))
{
	hashtab *self;

	DEBUG_ENTRY("newhash()");

	self = (hashtab *) MALLOC(sizeof(hashtab));
	self->freedata = freedata;
	self->nelem		 = 0;
	self->size		 = 0x1;
	self->fullmask = 0x1;
	self->frontmask = 0x0;
	self->space = 128;
	self->hashfunction = hashfunction;
	self->tab = (entry **) MALLOC(self->space*sizeof(entry *));
	self->tab[0] = NULL;

	return self;
}

void freehash(hashtab *table)
{
	entry *s, *t;
	unsigned long i;

	DEBUG_ENTRY("freehash()");

	for(i=0;i<table->size;i++) 
	{
		s = table->tab[i];
		while (s != NULL) 
		{
			t = s->next;
			if(table->freedata != NULL) 
			{
				table->freedata(s->data.p);
			}
			free(s);
			s = t;
		}
	}

	free(table->tab);
	free(table);
}

data_u *addentry(const void *key, size_t lkey, hashtab *table, int *exists)
{
	unsigned long hk, i;
	entry *s;
	void *p;

	DEBUG_ENTRY("addentry()");

	s = lookup_key(key,&lkey,table,&hk);
	if(s) 
	{ 
		*exists = 1;
		return &(s->data);
	}

	*exists = 0;
	p = MALLOC(sizeof(entry)+lkey);
	s = (entry *) p;
	s->hashval = hk;
	s->data.lkey = lkey;
	s->data.key = (void *)(((char *)p)+sizeof(entry));
	memcpy(s->data.key,key,lkey);

	i = getbin(hk,table);
	s->next = table->tab[i];
	table->tab[i] = s;

	table->nelem++;

	extend(table);

	return &(s->data);
}

data_u *lookup(const void *key, size_t lkey, const hashtab *table)
{
	unsigned long hk;
	entry *s;

	if(table->nelem == 0)
		return NULL;
	s = lookup_key(key,&lkey,table,&hk);
	if(s == NULL)
		return NULL;
	return &(s->data);
}

int maxchain(const hashtab *table)
{
	unsigned long i, l, max=0;
	entry *s;

	DEBUG_ENTRY("maxchain()");

	for(i=0;i<table->size;i++) 
	{
		l = 0;
		for(s = table->tab[i];s != NULL;s=s->next) 
		{
			l++;
		}
		if(l > max)
			max = l;
	}

	return max;
}

/* Internal utility functions */
/* Bernstein hash, see:
	 http://www.cse.yorku.ca/~oz/hash.html
	 http://www.burtleburtle.net/bob/hash/doobs.html 
*/
enum {HASHINIT=5381,HASHMUL=33}; /* Bernstein */
/* enum {HASHINIT=0,HASHMUL=131}; */
STATIC unsigned long hashfunction(const void *t, const size_t len)
{ 
	size_t i;
	unsigned long h = HASHINIT;
	unsigned char *s = (unsigned char *)t;

	DEBUG_ENTRY("hashfunction()");

	for( i=0; i<len; i++ ) 
		h = HASHMUL*h + s[i]; 

	return( h );
} 

STATIC entry *lookup_key(const void *key, size_t *lkey, const hashtab *table, 
			 unsigned long *hk)
{
	unsigned long i;
	entry *s;

	DEBUG_ENTRY("lookup_key()");

	if(*lkey == 0)
		*lkey = strlen((char *) key)+1;

	*hk = table->hashfunction(key,*lkey);
	i = getbin(*hk,table);

	/* Look through list within bin */
	for(s = table->tab[i];s!=NULL;s=s->next) 
	{
		/* Check for match -- full hash value will likely distinguish */
		if(*hk == s->hashval &&
				*lkey == s->data.lkey	 &&
				!memcmp(key,s->data.key,*lkey)) 
		{
			return s;
		}
	}
	return NULL;
}

STATIC void extend(hashtab *table)
{
	unsigned long move, last, i, j;
	entry *t, *s;

	DEBUG_ENTRY("extend()");
	last = table->size;
	table->size++;
	if(last > table->fullmask) 
	{	 /* Extend table when full */
		table->frontmask = table->fullmask;
		table->fullmask <<= 1;
		table->fullmask |= 1;
		if(table->fullmask+1 > table->space) 
		{
			table->space = table->fullmask+1;
			table->tab = (entry **) 
				REALLOC(table->tab,(table->space)*sizeof(entry *));
		}
	}

	/* Split next bin in front part */
	i = last & table->frontmask;	
	t = table->tab[i];
	table->tab[last] = table->tab[i] = NULL;
	move = table->frontmask ^ table->fullmask;
	while (t) 
	{ /* Go through list and sort between [last] and [i] */
		if(t->hashval & move) 
		{
			j = last;
		} 
		else 
		{
			j = i;
		}
		s = t->next;
		t->next = table->tab[j];
		table->tab[j] = t;
		t = s;
	}
}

STATIC unsigned long getbin(unsigned long hk, const hashtab *table)
{
	unsigned long i;

	DEBUG_ENTRY("getbin()");
	i = hk & table->fullmask; 
	if(i >= table->size)
		i &= table->frontmask;
	assert (i < table->size && i < table->space);
	return i;
}

/* Copy data from hash into list, 
	 optionally selected according to a masking function */
unsigned long makelist(const hashtab *table, data_u **list, 
											 const unsigned long nlist, int (*maskfun)(data_u *dat))
{
	unsigned long i, n;
	entry *s;

	DEBUG_ENTRY("makelist()");

	n = 0;
	for(i=0;i<table->size;i++) 
	{
		for(s = table->tab[i];s != NULL;s = s->next) 
		{
			if(maskfun == NULL || maskfun(&(s->data))) 
				list[n++] = &(s->data);
			if(n > nlist)
			{
				fprintf(ioQQQ,"Too many list items for array provided in makelist\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	return n;
}
unsigned long makeplist(const hashtab *table, void **list, 
											 const unsigned long nlist, int (*maskfun)(data_u *dat))
{
	unsigned long i, n;
	entry *s;

	DEBUG_ENTRY("makeplist()");
	n = 0;
	for(i=0;i<table->size;i++) 
	{
		for(s = table->tab[i];s != NULL;s = s->next) 
		{
			if(maskfun == NULL || maskfun(&(s->data))) 
				list[n++] = s->data.p;
			if(n > nlist)
			{
				fprintf(ioQQQ,"Too many list items for array provided in makeplist\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	return n;
}

#ifdef TEST
main()
{
	hashtab *table;
	data_u *data;
	int i, exists;

	char strings[][15] = {"Hello","Goodbye","Thirteen","One",
			"Two","Three","Four","Five","Six",
			"Seven","Eight"};

	table = newhash(NULL);
	for(i=0;i<sizeof(strings)/15;i++) 
	{
		data = addentry(strings[i],0,table,&exists);
		data->i = i;
	}

	for(i=0;i<sizeof(strings)/15;i++) 
	{ 
		data = lookup(strings[i],0,table);
		if(data) 
		{
			printf("Found data %d\n",data->i);
		}
		else 
		{
			printf("Couldn't find data\n");
		}
	}
	data = lookup("Banana",0,table);
	if(data) 
	{
		printf("Found data %d\n",data->i);
	} 
	else 
	{
		printf("Couldn't find data (as expected)\n");
	}
	printf("Longest chain is %d\n",maxchain(table));
}
#endif
