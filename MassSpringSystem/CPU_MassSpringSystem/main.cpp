
#include <iostream>
#include <vector>
#include <inttypes.h>
#include <memory>
#include <sstream>

#include <gl/freeglut.h>
#include "MassSpringMesh.h"

const char* WINDOW_TITLE = "SGPU - CPU - Mass-spring System";
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
GLint g_glutWindowIdentifier;

float FPS = 30.0f;
float INV_FPS = 1.0f / FPS;
const static unsigned int TIMERMSECS = static_cast<unsigned int>(INV_FPS * 1000);
bool paused = true;
std::shared_ptr<MassSpringMesh> mesh = nullptr;

void sinc(float x, float y, float& z) {
    z = std::sqrt(x*x + y*y);
}

void GenerateSincWave(float magnitude = 0.1f, float height = 1.0f) {
    std::vector<Particle>& nodes = mesh->getNodes();
    for ( std::size_t i = 0; i < nodes.size(); i++ ) {
        float x = nodes[i].position.getX();
        float z = nodes[i].position.getZ();

        float y;
        sinc(x, z, y);
        y = std::abs(y);
        y *= magnitude;
        y += height;
        mesh->setVelocity(i, Vector3f(0.0f, y, 0.0f));
    }

    mesh->setRestVelocitiesFromCurrent();
}

void g_init() {
	glClearColor(0.16f, 0.16f, 0.16f, 1.0f);

    std::cout << "[Control 1] Pause/Resume: Space Bar" << std::endl;
    std::cout << "[Control 2] Induce Wave: 's' key" << std::endl;

    mesh = std::make_shared<MassSpringMesh>();
    mesh->load("cloth.obj", 64.0f, 0.5f, Vector3f(0.0f, -9.81f, 0.0f));
    GenerateSincWave(1.0f, 1.0f);
}

void g_glutReshapeFunc(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, ((GLfloat)width/(GLfloat)height), 0.0, 1000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glutPostRedisplay();
}

void Project(float x, float y, float z, int& outx, int& outy) {
	double modelview[16];
	double projection[16];
	int viewport[4];
 
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

    double ox, oy, oz;
	gluProject(x, y, z, modelview, projection, viewport, &ox, &oy, &oz);
    outx = static_cast<int>(ox);
    outy = static_cast<int>(viewport[3] - oy);
}

void Draw2DText(const std::string& text, float x, float y, float z, float r, float g, float b) {
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

    int outx, outy;
    Project(x, y, z, outx, outy);
	glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		glOrtho(0.0, viewport[2], viewport[3], 0.0f, -1.0, 1.0); 
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadIdentity();
			glColor3f(r, g, b);
            glRasterPos2i(outx, outy);
			for ( int i = 0; i < text.length(); i++ )
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void DrawForces() {
    const std::vector<Particle>& nodes = mesh->getNodes();
    for ( std::size_t i = 0; i < nodes.size(); i++ ) {
        float x = nodes[i].position.getX();
        float y = nodes[i].position.getY();
        float z = nodes[i].position.getZ();

        std::stringstream stream;
        stream << "[" << nodes[i].force.getX() << " " << nodes[i].force.getY() << " " << nodes[i].force.getZ() << "]";
        Draw2DText(stream.str(), x, y, z, 1.0f, 1.0f, 1.0f);
    }
}

void DrawGrid() {
    glColor3f(0.18f, 0.18f, 0.18f);
    glBegin(GL_TRIANGLES);
        glVertex3f(-2.5f, 0.0f, -2.5f);
        glVertex3f(-2.5f, 0.0f, 2.5f);
        glVertex3f(2.5f, 0.0f, 2.5f);
        glVertex3f(-2.5f, 0.0f, -2.5f);
        glVertex3f(2.5f, 0.0f, 2.5f);
        glVertex3f(2.5f, 0.0f, -2.5f);
    glEnd();

    glColor3f(0.3f, 0.3f, 0.3f);
    glLineWidth(1.2f);
    glBegin(GL_LINES);
    for (GLfloat i = -2.5; i <= 2.5; i += 0.25) {
      glVertex3f(i, 0, 2.5); glVertex3f(i, 0, -2.5);
      glVertex3f(2.5, 0, i); glVertex3f(-2.5, 0, i);
    }
    glEnd();
}

bool RenderMassSpringMesh(const std::shared_ptr<MassSpringMesh>& mesh) {
    if ( mesh == nullptr ) return false;
    if ( mesh->getNodes().size() == 0 ) return false;
    
    const std::vector<Particle>& nodes = mesh->getNodes();
    const std::vector<Edge>& edges = mesh->getEdges();

    glPointSize(4.0f);
    glColor3f(0.2f, 0.6f, 1.0f);
    glBegin(GL_POINTS);
        for ( std::size_t i = 0; i < nodes.size(); i++ ) {
            glVertex3fv(nodes[i].position);
        }
    glEnd();

    glBegin(GL_LINES);
        for ( std::size_t i = 0; i < edges.size(); i++ ) {
            glVertex3fv(nodes[edges[i].p].position);
            glVertex3fv(nodes[edges[i].q].position);
        }
    glEnd();

    return true;
}

void g_glutDisplayFunc() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    glPushMatrix();
        gluLookAt(8.0, 3.0, 8.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

        DrawGrid();
        RenderMassSpringMesh(mesh);

    glPopMatrix();

	glFlush();
}

void g_glutKeyboardFunc(unsigned char key, int x, int y) {
    if ( key == ' ' ) paused = !paused;
    if ( key == 's' ) GenerateSincWave(1.0f, 1.0f);
}

void update(int value) {
    glutTimerFunc(TIMERMSECS, update, 0);
    if ( paused == false ) {
        mesh->update(INV_FPS);
        mesh->resolveCollisions();
    }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	g_glutWindowIdentifier = glutCreateWindow(WINDOW_TITLE);

	glutDisplayFunc(g_glutDisplayFunc);
	glutReshapeFunc(g_glutReshapeFunc);
    glutKeyboardFunc(g_glutKeyboardFunc);

    glutTimerFunc(TIMERMSECS, update, 0);

	g_init();

	glutMainLoop();
	return 0;
}
