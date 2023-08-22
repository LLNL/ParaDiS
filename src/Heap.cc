/**************************************************************************
 *
 *      Module:       Heap.c
 *      Description:  Contains functions for maintaining an integer "heap".
 *                    The heap is a binary tree. The rule for a heap is
 *                    that the value of a cell is greater than or equal
 *                    to the value of its parent cell. This applies
 *                    recursively throughout the tree (heap). This ensures
 *                    that the smallest number is in the root cell. 
 *                    This tree is stored in a linear array such that
 *                    the children of any element k are located at
 *                    elements 2k+1 and 2k+2.
 *
 *                    To add an element, one can put it at the end of
 *                    the array (which would mean adding a new node at
 *                    the bottom of the tree). At this point the heap rule
 *                    may not be fulfilled. If one element was changed
 *                    or added, it is only around that element that the rule
 *                    can be broken. So we can compare the value of the new
 *                    node to that if its parent. If the rule is not
 *                    fulfilled, we swap them. If we made a swap, we
 *                    recursively apply the same test to the parent node and
 *                    possible switch. If no swap is made, the heap condition
 *                    is fulfilled everywhere and we are done. If we reach
 *                    the root node, we are done, since the root node has no
 *                    parent, and the heap condition is fullfilled everywhere
 *                    below. This process is called bubble up, since an added
 *                    small value will wiggle its way upward like a bubble in
 *                    water.
 *
 *                    If the top element is removed, we have an empty spot
 *                    in the array that needs to be filled. The simplest way
 *                    is to take the last element and put it in the empty
 *                    first spot. Then a bubble down is performed. Here a
 *                    cell value is compared to its children values, and if
 *                    a childs value is smaller, a swap is made, and the
 *                    process continues recursively downward
 *
 *      Functions:
 *          BubbleDown()
 *          BubbleUp()
 *          HeapAdd()
 *          HeapInit()
 *          HeapPeek()
 *          HeapRemove()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

/*
 *      Define the number of elements by which to increase the heap
 *      size when more space is needed.
 */
#define HEAP_SIZE_INCR 50


/*------------------------------------------------------------------------
 *
 *      Function:     BubbleUp
 *      Description:  Shift the value at the specified location in
 *                    the heap up the heap to its proper position
 *
 *      Arguments:
 *          heap      Pointer to the array containing the heap
 *          index     Index in the array <heap> of the data element
 *                    that needs to be shifted upward in the heap
 *
 *----------------------------------------------------------------------*/
static void BubbleUp(int *heap, int index) 
{
        int k, last, tmp;

        if (index == 0) return;

        last = index;

        while (1) {
            k = (last - 1) >> 1;
            if (k >= 0) {
                if (heap[last] < heap[k]) {
                    tmp = heap[k];
                    heap[k] = heap[last];
                    heap[last] = tmp;
                    last = k;
                    continue;
                }
            }
            break;
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     BubbleDown
 *      Description:  Shift the value at the specified location in
 *                    the heap down the heap to its proper position
 *
 *      Arguments:
 *          heap      Pointer to the array containing the heap
 *          heapCnt   Number of heap elements currently in use
 *          index     Index in the array <heap> of the data element
 *                    that needs to be shifted downward in the heap
 * 
 *----------------------------------------------------------------------*/
static void BubbleDown(int *heap, int heapCnt, int index)
{
        int k, left, right, tmp;

        k = index;

        while (1) {
            left = k * 2 + 1;
            right = left + 1;
            if (left < heapCnt) {
                if (right < heapCnt) {
                    left = (heap[right] < heap[left] ? right : left);
                }
                if (heap[left] < heap[k]) {
                    tmp = heap[k];
                    heap[k] = heap[left];
                    heap[left] = tmp;
                    k = left;
                    continue;
                }
            }
            break;
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     HeapAdd
 *      Description:  Insert the provided value into the heap preserving
 *                    the condition that all elements in the tree are
 *                    of value greater or equal to their parent elements.
 *                    If the heap is out of space, the heap size will
 *                    be increased and the new heap pointer returned to
 *                    the caller.
 *
 *      Arguments:
 *          heap      Pointer to the heap array pointer.  If the heap
 *                    needs to grow in order to accomodate the new
 *                    value, the heap will be reallocated and the 
 *                    new pointer returned to the caller in this value.
 *          heapSize  Pointer to the size of the currently allocated
 *                    heap (as a count of integers).  If the heap needs
 *                    to grow in order to accomodate the new value, the
 *                    new heapSize will be returned to the caller in this
 *                    location.
 *          heapCnt   Pointer to the number of elements of the heap 
 *                    currently in use.  The value at this location
 *                    will be incremented to account for the newly
 *                    added element.
 *          value     Integer value to be added to the heap.
 * 
 *----------------------------------------------------------------------*/
void HeapAdd(int **heap, int *heapSize, int *heapCnt, int value)
{
        int newCnt;

        newCnt = *heapCnt + 1;

/*
 *      If we've run out of space in the heap, increase the heap size;
 */
        if (newCnt >= *heapSize) {
            *heapSize += HEAP_SIZE_INCR;
            *heap = (int *)realloc(*heap, *heapSize * sizeof(int));
        }

        (*heap)[*heapCnt] = value;
        BubbleUp(*heap, *heapCnt);

        *heapCnt = newCnt;

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     HeapRemove
 *      Description:  Remove the smallest value from the specified heap
 *                    and return it to the caller.
 *
 *      Arguments:
 *          heap      Pointer to the heap array.
 *          heapCnt   Pointer to the number of elements of the heap 
 *                    currently in use.  The value at this location
 *                    will be decremented to account for the newly
 *                    removed element.
 * 
 *      Returns:  lowest available recycled tag index on success, or
 *                -1 if no recylced tags available.
 *
 *----------------------------------------------------------------------*/
int HeapRemove(int *heap, int *heapCnt)
{
        int value;

/*
 *      If the heap in not empty, return the first (smallest) element,
 *      otherwise, return a -1.
 */
        if (*heapCnt > 0) {
            value = heap[0];
            *heapCnt -= 1;
            heap[0] = heap[*heapCnt];
            BubbleDown(heap, *heapCnt, 0);
        } else {
            value = -1;
        }

        return(value);
}


/*------------------------------------------------------------------------
 *
 *      Function:     HeapPeek
 *      Description:  Take a peek at the smallest available value on the
 *                    heap and return it to the caller without removing it
 *                    from the heap.
 *
 *      Arguments:
 *          heap      Pointer to the heap array.
 *          heapCnt   Number of elements of the heap currently in use.
 * 
 *      Returns:  lowest available value on the heap on success, or
 *                -1 if the heap is empty.
 *
 *      NOTE: Might be better just to define a macro for this rather
 *            than a function.
 *
 *----------------------------------------------------------------------*/
int HeapPeek(int *heap, int heapCnt)
{
    return((heapCnt > 0) ? heap[0] : -1);
}
