#ifndef GENERIC_LIST_H
#define GENERIC_LIST_H

#define INITIAL_SIZE 100000

#include <cstdlib>
#include <stdexcept>
#include <string.h>

template <class T>

class generic_list {

public:
	generic_list(void);

	~generic_list(void);

	int append(T new_elem);
	int remove(int target_elem);
	T get_elem(int index);
	int set_elem(int index, T new_elem);
	int get_num_elems();

private:
	int numElems;
	int maxElems;
	T *theElems;
};
template <class T>
generic_list<T>::generic_list()
{
    numElems=0;
    maxElems = INITIAL_SIZE;
    theElems = (T *)malloc(maxElems * sizeof(T));
}

template <class T>
generic_list<T>::~generic_list(void) 
{
    free(theElems);
}

template <class T>
int generic_list<T>::set_elem(int index, T new_elem)
{
    if (index > numElems-1)
	return -1;

    // insert the cell
    theElems[index]=new_elem;

    return 0;
}
template <class T>
int generic_list<T>::append(T new_elem)
{
    if (numElems + 1 > maxElems)
    {
	// allocate double current storage
	T *tmpList = (T *)malloc (maxElems*2*sizeof(T)); 

	// copy existing data
	memcpy(tmpList, theElems, maxElems*sizeof(T));

	//fprintf(stderr,"Realloc of list %x (new max=%d)\n",
	 //   (void *)this,maxElems*2);

	// free old data
	free(theElems);

	// switch to new list storage
	theElems=tmpList;

	maxElems=maxElems*2;
    }
    // insert the cell
    theElems[numElems]=new_elem;
    numElems++;

    return 0;
}

template <class T>
int generic_list<T>::remove(int target_elem)
{
    if (target_elem >= 0 && target_elem < numElems)
    {
	// reclaim the space by moving all subsequent cell points down
	for (int i=target_elem; i < numElems; i++)
	    theElems[i]=theElems[i+1];

	numElems--;
	return 0;
    }
    return -1;
}

template <class T>
T generic_list<T>::get_elem(int index)
{
    if (index > numElems-1)
	return NULL;

    return theElems[index];
}

template <class T>
int generic_list<T>::get_num_elems()
{
    return numElems;
}
#endif
