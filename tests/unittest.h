#pragma once

#include "qdp.h"
#include <vector>

#include "cli_args.h"

#include <qphix/qphix_config.h>
using namespace QDP;
using namespace std;

namespace Assertions
{
template <typename T>
inline void assertEquals(const T &t1, const T &t2)
{
  if (t1 != t2) {
    throw std::exception();
  }
}

inline void fail(const std::string &message)
{
  QDPIO::cout << "\t " << message << "  ";
  throw std::exception();
}

template <typename T>
inline void assertNotEquals(const T &t1, const T &t2)
{
  if (t1 == t2) {
    throw std::exception();
  }
}

inline void assertion(bool b)
{
  if (!b) {
    throw std::exception();
  }
}

// Add other assertions here
};

// A test case - provides a run() method and a getName() to identify itself
// Strictly speaking the getName() is probably not necessary
class TestCase
{
 private:
 public:
  virtual void run(void) = 0;
  virtual ~TestCase() {}
};

// A test fixture - that does extra set up before and after the test
class TestFixture : public TestCase
{
 public:
  virtual ~TestFixture() {}
  virtual void setUp() {}
  virtual void tearDown() {}
  virtual void runTest() {}

  void run(void)
  {
    setUp();
    runTest();
    tearDown();
  }
};

// A runner class, which holds a bunch of tests
// Essentially this is a testCase,  but it is special
// As it sets up QDP as well -- can probably fudge things
// to get the lattice dimensions from argc, and argv
class TestRunner : public TestCase
{
 private:
  int num_success;
  int num_failed;
  int num_unexpected_failed;
  int num_tried;

  enum Success { SUCCESS, FAIL, ERROR };

  struct TestItem {
    TestItem(TestCase *t, const std::string &n) : test(t), name(n){};
    ~TestItem() { delete test; }
    TestCase *test;
    const std::string name;
    Success success;
  };

  std::vector<TestItem *> tests;

  CliArgs args_;

 public:
  TestRunner(int *argc, char ***argv)
      : num_success(0), num_failed(0), num_unexpected_failed(0), num_tried(0)
  {
    QDP_initialize(argc, argv);

    args_ = processArgs(*argc, *argv);

    multi1d<int> nrow(Nd);
    for (int i = 0; i < Nd; ++i) {
      nrow[i] = args_.nrow_in[i];
    }
    Layout::setLattSize(nrow);
    Layout::create();

    omp_set_num_threads( (args_.NCores + args_.NCommCores) * args_.Sy * args_.Sz);
  }

  CliArgs &args() { return args_; }

  void run(void)
  {
    for (unsigned i = 0; i != tests.size(); i++) {
      run(*(tests[i]));
    }
  }

  // Extra Features: Add a test
  void addTest(TestCase *t, const std::string &name)
  {
    if (t != 0x0) {
      TestItem *ti = new TestItem(t, name);

      tests.push_back(ti);
    }
  }

  // Run a particular test
  void run(TestItem &t)
  {

    try {
      QDPIO::cout << "Running Test: " << t.name;
      num_tried++;
      t.test->run();
      t.success = SUCCESS;
    } catch (std::exception const &e) {
      QDPIO::cout << e.what() << std::endl;
      t.success = FAIL;
    } catch (...) {
      t.success = ERROR;
    }

    bool test_success = trueEverywhere(t.success == SUCCESS);
    if (test_success == true) {
      QDPIO::cout << "  OK" << std::endl;
    } else {

      if (trueOnNodes(t.success == ERROR) > 0) {
        QDPIO::cout << "  ERROR" << std::endl;
      } else {
        QDPIO::cout << "  FAIL" << std::endl;
      }
    }
  }

  // Compute the number of nodes on which condition is true
  int trueOnNodes(bool condition)
  {
    int summand;
    if (condition == true) {
      summand = 1;
    } else {
      summand = 0;
    }
    QDPInternal::globalSum(summand);
    return summand;
  }

  bool trueEverywhere(bool condition)
  {
    return (trueOnNodes(condition) == Layout::numNodes());
  }

  void summary()
  {
    QDPIO::cout << "Summary: " << num_tried << " Tests Tried" << std::endl;
    int success = 0;
    int failure = 0;
    int odd = 0;
    for (unsigned int i = 0; i < tests.size(); i++) {
      if (trueEverywhere(tests[i]->success == SUCCESS)) {
        success++;
      } else {
        failure++;
        if (trueOnNodes(tests[i]->success == ERROR) > 0) {
          odd++;
        }
      }
    }
    QDPIO::cout << "         " << success << " Tests Succeeded " << std::endl;
    QDPIO::cout << "         " << failure << " Tests Failed on some nodes"
                << std::endl;
    QDPIO::cout << "of which " << odd
                << " Tests Failed in Unexpected Ways on some nodes" << std::endl;
    if (failure > 0) {
      exit(1);
    }
  }

  ~TestRunner()
  {

    for (std::vector<TestItem *>::iterator i = tests.begin(); i != tests.end();
         i++) {
      delete (*i);
    }
    QDP_finalize();
  }

 private:
};
