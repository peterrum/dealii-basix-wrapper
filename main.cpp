#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <basix/finite-element.h>

#include <memory>

using namespace dealii;

template <int dim, int spacedim = dim>
class BasixFE : public FiniteElement<dim, spacedim>
{
public:
  BasixFE()
    : FiniteElement<dim, spacedim>(
        FiniteElementData<dim>(std::vector<unsigned int>() /*TODO*/,
                               ReferenceCells::get_simplex<dim>() /*TODO*/,
                               1 /*TODO*/,
                               1 /*TODO*/,
                               FiniteElementData<dim>::H1),
        std::vector<bool>(),
        std::vector<ComponentMask>())
  {}

  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override
  {
    return std::make_unique<BasixFE<dim, spacedim>>(*this);
  }

  std::string
  get_name() const override
  {
    std::ostringstream namebuf;
    namebuf << "BasixFE<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }

  UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override
  {
    Assert(false, ExcNotImplemented());

    (void)update_flags;

    return {};
  }

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    Assert(false, ExcNotImplemented());

    (void)update_flags;
    (void)mapping;
    (void)quadrature;
    (void)output_data;

    return {};
  }

  void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const
  {
    Assert(false, ExcNotImplemented());

    (void)cell;
    (void)cell_similarity;
    (void)quadrature;
    (void)mapping;
    (void)mapping_internal;
    (void)mapping_data;
    (void)fe_internal;
    (void)output_data;
  }

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    Assert(false, ExcNotImplemented());

    (void)cell;
    (void)face_no;
    (void)sub_no;
    (void)quadrature;
    (void)mapping;
    (void)mapping_internal;
    (void)mapping_data;
    (void)fe_internal;
    (void)output_data;
  }
};

/**
 * /home/munch/sw-basix/cmake-3.22.0-install/bin/cmake ..
 * -DBasix_DIR=/home/munch/sw-basix/basix-install/lib64/cmake/basix/
 * -Dxtl_DIR=/home/munch/sw-basix/basix/_deps/xtl-build
 * -Dxtensor_DIR=/home/munch/sw-basix/basix/_deps/xtensor-build
 * -DDEAL_II_DIR=../dealii-build
 */
int
main()
{
  auto _element = std::make_unique<basix::FiniteElement>(
    basix::create_element(basix::element::family::P,
                          basix::cell::type::triangle,
                          2,
                          basix::element::lagrange_variant::equispaced,
                          false));

  std::cout << _element->dim() << std::endl;

  for (int j = 0; j < _element->dim(); ++j)
    {
      for (unsigned int i = 0; i < 2; ++i)
        std::cout << _element->points()(j, i) << " ";
      std::cout << std::endl;
    }

  const auto &trafo = _element->base_transformations();

  std::cout << trafo.shape()[0] << std::endl;
  std::cout << trafo.shape()[1] << std::endl;
  std::cout << trafo.shape()[2] << std::endl;

  for (unsigned int i = 0; i < trafo.shape()[0]; ++i)
    {
      for (unsigned int j = 0; j < trafo.shape()[1]; ++j)
        {
          for (unsigned int k = 0; k < trafo.shape()[2]; ++k)
            std::cout << trafo(i, j, k) << " ";
          std::cout << std::endl;
        }
      std::cout << std::endl;
    }

  {
    const unsigned int dim = 2;

    Triangulation<dim> tria;

    GridGenerator::subdivided_hyper_rectangle_with_simplices(tria,
                                                             {1, 1},
                                                             {0.0, 0.0},
                                                             {1.0, 1.0});

    DoFHandler<dim> dof_handler(tria);

    BasixFE<dim> fe;

    dof_handler.distribute_dofs(fe);
  }

  return 0;
}
