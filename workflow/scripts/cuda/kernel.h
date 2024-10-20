#ifndef KERNEL_H
#define KERNEL_H


#include <stdio.h>
#include <iostream>
#include <math.h>


class Classifier {

private:
    int id;

public:
    Classifier(int id);
    int getId();
    void printId();

};

#endif

