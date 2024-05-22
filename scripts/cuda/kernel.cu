#include "kernel.h"

Classifier::Classifier(int classifierId) {
    id = classifierId;
}


int Classifier::getId() {
    return id;
}
void Classifier::printId(){
    printf("Classifier ID : %d\n", id);
}

int main() {

    Classifier classifier(46);
    classifier.printId();

    return 0;

}