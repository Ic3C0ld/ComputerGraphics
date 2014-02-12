#include "main.h"


int main(int argc, char* argv[])
{
	// initialize GLUT library state
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );

	glutInitWindowSize(1200, 800);
	glutInitWindowPosition(100, 50);

	// Create and label the main window
	glutCreateWindow("Course5");

	// Configure various properties of the OpenGL rendering context
	Setup();

	// The rendering function 
	glutDisplayFunc(Render);
	glutReshapeFunc(Resize);
	glutIdleFunc(Idle);

	glutKeyboardFunc(keyboardDown);
	glutKeyboardUpFunc(keyboardUp);
	glutSpecialFunc(specialKeyFunc);
	glutMouseFunc(mouseFunc);
	glutMotionFunc(mouseMotion);
	glutPassiveMotionFunc(mousePassiveMotion);

	//Enter main event handling loop
	glutMainLoop();
	return 0;
}