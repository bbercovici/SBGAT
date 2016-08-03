#include "utilities.h"
#include <iostream>
void emulate_key_press(const int keyCode) {
    CGEventRef eventRef = CGEventCreateKeyboardEvent (NULL, keyCode, true);
    CGEventPost(kCGSessionEventTap, eventRef);
    CFRelease(eventRef);
    std::cout << "r pressed" << std::endl;
}